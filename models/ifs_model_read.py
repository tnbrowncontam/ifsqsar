"""ifs_model_read.py
Reads model files from the development code for group contribution QSARs (IFS) developed by Trevor N. Brown.
Applies the QSARs and checks domain of applicability for single structures as SMILES or for a file of structures.
"""

try:
    # ob version 3.x
    from openbabel import openbabel as ob
except ImportError:
    # ob version 2.x
    import openbabel as ob
import numpy as np
import importlib


def _normcdfapprox(x):
    """Returns normal distribution CDF."""
    # approximation from:
    # Hector Vazquez-Leal, Roberto Castaneda-Sheissa, Uriel Filobello-Nino,
    # Arturo Sarmiento-Reyes, and Jesus Sanchez Orea, "High Accurate Simple
    # Approximation of Normal Distribution Integral," Mathematical Problems
    # in Engineering, vol. 2012, Article ID 124029, 22 pages, 2012.
    # doi:10.1155/2012/124029
    return 1./(np.exp(-x*358./23.+111*np.arctan(x*37./294.))+1.)


def _calculate_fragment_similarity(counts_array_i, counts_array_j, stdev_array, simil_cut=None):
    """Calculate Tanimoto similarity coefficient between two arrays of fragment counts."""
    # all b-values = 1 meaning:
    #  1) a-values * b-values = a-values
    #  2) b-values**2 = b-values
    #  3) sum(b-values) = # non-zero elements
    #  so only the number of non-zero elements needs to be known
    b = float(np.logical_or(counts_array_i, counts_array_j).sum())
    if b == 0.:
        tanimoto = 0.
    else:
        # find non-zero a-values
        amask1 = np.logical_and(counts_array_i, counts_array_j)
        a1 = amask1.sum()
        # at this point sum(a)/sum(b) is the maximum possible CSS for this
        # pair of chemicals, so if it's less than the similarity cutoff
        # then skip the rest of the calculation
        if simil_cut is not None and a1/b < simil_cut:
            tanimoto = 0.
        else:
            # some a-values will be equal to one because the fragment
            # counts are the same (a1s). only the number needs to be
            # known for similar reasons as the b-values
            amask2 = np.logical_and(amask1, counts_array_i != counts_array_j)
            a1 -= amask2.sum()
            # calculate non-zero and non-one a-values if any.
            # for fragments with stdev = 0, the a value is = 0,
            # so do not include these in the array
            zmask = stdev_array != 0
            aa = np.abs(counts_array_i[amask2*zmask] - counts_array_j[amask2*zmask]) / stdev_array[amask2*zmask]
            a = 1 - (2 * _normcdfapprox(aa) - 1)
            # calculate tanimoto => sum(ab) / [sum(a2) + sum(b2) - sum(ab)]
            asum = a1 + a.sum()
            tanimoto = asum / ((a1+(a**2).sum()) + b - asum)
    return tanimoto


class QSARModel:
    """Class that loads an arbitrary Model from a Model file and applies it
    to any molecules passed to it as an openbabel mol"""
    
    def __init__(self, model_module, model_name):
        """Save model name and create model namespace."""
        
        # define Model namespace
        self.model_namespace = None
        self.model_module = model_module
        self.model_name = model_name

    def load(self):
        """Import the model module."""
        self.model_namespace = importlib.import_module(self.model_module)
        self.model_namespace.smartslist = []
        for smarts in self.model_namespace.fragmentlist['smarts']:
            if smarts == b'intercept':
                self.model_namespace.smartslist.append('intercept')
            else:
                pattern = ob.OBSmartsPattern()
                pattern.Init(smarts.decode('utf-8'))
                self.model_namespace.smartslist.append(pattern)
        self.model_namespace.coefficientarray = np.mean(self.model_namespace.coefficientarrays, axis=1)
        if self.model_namespace.domain:
            if self.model_namespace.intercept:
                self.model_namespace.xtxi = np.linalg.inv(np.matmul(self.model_namespace.train_counts[:, 1:].T, self.model_namespace.train_counts[:, 1:]))
            else:
                self.model_namespace.xtxi = np.linalg.inv(np.matmul(self.model_namespace.train_counts.T, self.model_namespace.train_counts))
            self.model_namespace.neg_dom_check_init = []
            for s in range(self.model_namespace.neg_dom_check.shape[0]):
                smarts1, smarts2, description = self.model_namespace.neg_dom_check[s]
                pattern1 = ob.OBSmartsPattern()
                pattern1.Init(smarts1.decode('utf-8'))
                pattern2 = ob.OBSmartsPattern()
                pattern2.Init(smarts2.decode('utf-8'))
                self.model_namespace.neg_dom_check_init.append((pattern1, pattern2, description))

    def apply_model(self, molecule):
        """Take an openbabel molecule, apply the QSARModel and return the result."""
        # check if model has been loaded
        if self.model_namespace is None:
            self.load()
        # add or delete hydrogens depending on model
        if self.model_namespace.molecule_format == 'old_format':
            molecule.AddHydrogens()
        else:
            molecule.DeleteHydrogens()
        # get fragment counts
        fragment_counts = []
        for smarts in self.model_namespace.smartslist:
            if smarts == 'intercept':
                fragment_counts.append(1)
            else:
                smarts.Match(molecule)
                fragment_counts.append(len(smarts.GetUMapList()))
        fragment_counts = np.array(fragment_counts)
        # apply multiple linear regression qsar
        if self.model_namespace.model_type == 'MLR':
            # apply qsar
            prediction = (fragment_counts * self.model_namespace.coefficientarray).sum()
            error = np.nan
            if self.model_namespace.domain:
                # calculate CSS
                topn = 5
                topgroup = []
                mintop = 0.
                if self.model_namespace.intercept:
                    ibegin = 1
                else:
                    ibegin = 0
                for t in range(self.model_namespace.train_counts.shape[0]):
                    fragsim = _calculate_fragment_similarity(fragment_counts[ibegin:],
                                                             self.model_namespace.train_counts[t, ibegin:],
                                                             self.model_namespace.fragmentlist['fragstdev'][ibegin:],
                                                             mintop)
                    topgroup.append((fragsim, 1-self.model_namespace.datalist['value_similarity'][t]))
                    topgroup.sort(reverse=True)
                    if len(topgroup) > topn:
                        for i in reversed(range(topn, len(topgroup))):
                            if topgroup[i][0] != topgroup[topn-1][0]:
                                topgroup.pop(i)
                    mintop = topgroup[-1][0]
                css = 1
                for i in range(topn):
                    css *= (topgroup[i][0]*(1-topgroup[i][1]))**0.5
                css = css**(1/float(topn))
                # calculate leverage
                if self.model_namespace.intercept:
                    leverage = np.matmul(np.matmul(fragment_counts[1:], self.model_namespace.xtxi), fragment_counts[1:].T)
                else:
                    leverage = np.matmul(np.matmul(fragment_counts, self.model_namespace.xtxi), fragment_counts.T)
                # set uncertainty level and note
                ul = 0
                error = self.model_namespace.warn_0_error
                note = []
                if self.model_namespace.intercept:
                    sumcounts = fragment_counts[1:].sum()
                else:
                    sumcounts = fragment_counts.sum()
                if sumcounts == 0:
                    ul = 4
                    error = self.model_namespace.warn_4_error
                    note.append('no fragment overlap with training dataset')
                elif leverage >= 1:
                    ul = 3
                    error = self.model_namespace.warn_3_error
                    note.append('leverage > 1')
                elif css <= self.model_namespace.css_cutoff_1 or leverage >= self.model_namespace.leverage_cutoff_1:
                    ul = 2
                    error = self.model_namespace.warn_2_error
                    if css < self.model_namespace.css_cutoff_1:
                        note.append('out of domain')
                    if leverage > self.model_namespace.leverage_cutoff_1:
                        note.append('structural outlier')
                elif css <= self.model_namespace.css_cutoff_0 or leverage >= self.model_namespace.leverage_cutoff_0:
                    ul = 1
                    error = self.model_namespace.warn_1_error
                    if css < self.model_namespace.css_cutoff_0:
                        note.append('low similarity')
                    if leverage > self.model_namespace.leverage_cutoff_0:
                        note.append('high leverage')
                # negative domain check for atom type violations
                violations = []
                for pattern1, pattern2, description in self.model_namespace.neg_dom_check_init:
                    pattern1.Match(molecule)
                    pattern2.Match(molecule)
                    if len(pattern1.GetUMapList()) != len(pattern2.GetUMapList()):
                        violations.append(description)
                if len(violations) > 0:
                    ul = 5
                    error = self.model_namespace.warn_5_error
                    for v in violations:
                        note.append(v)
                # check bounds
                if self.model_namespace.lower_bound and prediction < self.model_namespace.min_train:
                    ul = 6
                    note.append('prediction (' + str(round(prediction, self.model_namespace.round_digits)) + ') less than smallest value in training set')
                    prediction = self.model_namespace.min_train
                if self.model_namespace.upper_bound and prediction > self.model_namespace.max_train:
                    ul = 6
                    note.append('prediction (' + str(round(prediction, self.model_namespace.round_digits)) + ') greater than largest value in training set')
                    prediction = self.model_namespace.max_train
                post_proc_prediction, post_proc_error = self.model_namespace.post_processing(prediction, error)
                return post_proc_prediction, ul, post_proc_error, ', '.join(note)
            else:
                post_proc_prediction, post_proc_error = self.model_namespace.post_processing(prediction, error)
                return post_proc_prediction, np.nan, post_proc_error, ''

