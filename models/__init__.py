"""
ifsqsar/models subpackage
developed by Trevor N. Brown
Stores all data and model-specific code for QSARs as python modules and implements a generic API for accessing them
"""

from openbabel import openbabel as ob
import numpy as np
import importlib


def _normcdfapprox(x):
    """Returns normal distribution CDF"""
    # approximation from:
    # Hector Vazquez-Leal, Roberto Castaneda-Sheissa, Uriel Filobello-Nino,
    # Arturo Sarmiento-Reyes, and Jesus Sanchez Orea, "High Accurate Simple
    # Approximation of Normal Distribution Integral," Mathematical Problems
    # in Engineering, vol. 2012, Article ID 124029, 22 pages, 2012.
    # doi:10.1155/2012/124029
    return 1. / (np.exp(-x * 358. / 23. + 111 * np.arctan(x * 37. / 294.)) + 1.)


def _calculate_fragment_similarity(counts_array_i, counts_array_j, stdev_array, simil_cut=None):
    """Calculate Tanimoto similarity coefficient between two arrays of fragment counts"""
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
        if simil_cut is not None and a1 / b < simil_cut:
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
            aa = np.abs(counts_array_i[amask2 * zmask] - counts_array_j[amask2 * zmask]) / stdev_array[amask2 * zmask]
            a = 1 - (2 * _normcdfapprox(aa) - 1)
            # calculate tanimoto => sum(ab) / [sum(a2) + sum(b2) - sum(ab)]
            asum = a1 + a.sum()
            tanimoto = asum / ((a1 + (a ** 2).sum()) + b - asum)
    return tanimoto


class QSARModel:
    """Class that loads a QSAR stored as a python module and applies it
    to molecules passed to it as IFSMols (subclass of openbabel mols)"""

    def __init__(self, model_module, model_name, version):
        """Save model name and create model namespace"""
        # define Model namespace
        self.model_namespace = None
        self.ismixture = False
        self.model_module = model_module
        self.model_name = model_name
        self.version = version
        self.super_models = []
        self.default_stored = {}

    def __str__(self):
        return self.model_name

    def __repr__(self):
        return self.__str__()

    def load(self):
        """Import the QSAR python module"""
        if self.model_namespace is not None:
            return
        self.model_namespace = importlib.import_module(self.model_module)
        self.model_namespace.smartslist = []
        for smarts in self.model_namespace.fragmentlist['smarts']:
            if smarts == b'intercept':
                self.model_namespace.smartslist.append('intercept')
            elif smarts == b'sssr':
                self.model_namespace.smartslist.append('sssr')
            elif smarts == b'MW':
                self.model_namespace.smartslist.append('MW')
            else:
                pattern = ob.OBSmartsPattern()
                pattern.Init(smarts.decode('utf-8'))
                self.model_namespace.smartslist.append(pattern)
        self.model_namespace.coefficientarray = np.mean(self.model_namespace.coefficientarrays, axis=1)
        if self.model_namespace.domain:
            if self.model_namespace.intercept:
                self.model_namespace.xtxi = np.linalg.inv(
                    np.matmul(self.model_namespace.train_counts[:, 1:].T, self.model_namespace.train_counts[:, 1:]))
            else:
                self.model_namespace.xtxi = np.linalg.inv(
                    np.matmul(self.model_namespace.train_counts.T, self.model_namespace.train_counts))
            self.model_namespace.neg_dom_check_init = []
            for s in range(self.model_namespace.neg_dom_check.shape[0]):
                smarts1, smarts2, description = self.model_namespace.neg_dom_check[s]
                pattern1 = ob.OBSmartsPattern()
                pattern1.Init(smarts1.decode('utf-8'))
                pattern2 = ob.OBSmartsPattern()
                pattern2.Init(smarts2.decode('utf-8'))
                self.model_namespace.neg_dom_check_init.append((pattern1, pattern2, description.decode('utf-8')))
        # backup the stored data for reset and restore
        self.default_stored = self.model_namespace.stored.copy()

    def set_stored(self, normsmiles, value, ul='U', error=np.nan, ulnote='user value', citation='user value', units=None, endpoint=None):
        """Add a stored value with user-entered data"""
        # check if model has been loaded
        if self.model_namespace is None:
            self.load()
        # check that UL is one of the permitted values
        if type(ul) is not int and type(ul) is not str:
            raise ValueError('User-defined UL must be E, U, int 0-6 or combination thereof')
        elif type(ul) is int and ul not in [0, 1, 2, 3, 4, 5, 6]:
            raise ValueError('User-defined UL must be E, U, int 0-6 or combination thereof')
        elif type(ul) is str:
            counts = []
            for c in ['E', 'U', '0', '1', '2', '3', '4', '5', '6']:
                counts.append(ul.count(c))
                if counts[-1] > 1:
                    raise ValueError('User-defined UL must be E, U, int 0-6 or combination thereof')
            if len(ul) - sum(counts) > 0:
                raise ValueError('User-defined UL must be E, U, int 0-6 or combination thereof')
        # set units and endpoint from model namespace if not provided with input
        if units is None:
            units = self.model_namespace.units
        if endpoint is None:
            endpoint = self.model_namespace.endpoint
        # add user value to stored data
        self.model_namespace.stored[normsmiles] = (value, ul, error, ulnote, citation, units, endpoint)

    def load_stored(self, normsmiles, propagatedown=False, propagateup=False):
        """Remove a specific stored value"""
        # check if model has been loaded
        if self.model_namespace is None:
            self.load()
        if normsmiles in self.default_stored:
            self.model_namespace.stored[normsmiles] = self.default_stored[normsmiles]
        if propagateup:
            for qsar in self.super_models:
                qsar.load_stored(normsmiles, propagateup=propagateup)

    def remove_stored(self, normsmiles, propagatedown=False, propagateup=False):
        """Remove a specific stored value"""
        # check if model has been loaded
        if self.model_namespace is None:
            self.load()
        if normsmiles in self.model_namespace.stored:
            self.model_namespace.stored.pop(normsmiles)
        if propagateup:
            for qsar in self.super_models:
                qsar.remove_stored(normsmiles, propagateup=propagateup)

    def reset_stored(self, propagatedown=False, propagateup=False):
        """Reset the stored data to default, erasing all predicted and user data"""
        # check if model has been loaded
        if self.model_namespace is None:
            self.load()
        self.model_namespace.stored = self.default_stored.copy()
        if propagateup:
            for qsar in self.super_models:
                qsar.reset_stored(propagateup=propagateup)

    def restore_stored(self, propagatedown=False, propagateup=False):
        """Reset the stored data to default, preserving predicted and user data not in the default data"""
        # check if model has been loaded
        if self.model_namespace is None:
            self.load()
        self.model_namespace.stored.update(self.default_stored)
        if propagateup:
            for qsar in self.super_models:
                qsar.restore_stored(propagateup=propagateup)

    def erase_user_stored(self, propagatedown=False, propagateup=False):
        """Erase all stored data flagged as user ('U')"""
        # check if model has been loaded
        if self.model_namespace is None:
            self.load()
        for key in list(self.model_namespace.stored.keys()):
            if self.model_namespace.stored[key][1] == 'U':
                self.model_namespace.stored.pop(key)
        if propagateup:
            for qsar in self.super_models:
                qsar.erase_user_stored(propagateup=propagateup)

    def erase_experimental_stored(self, propagatedown=False, propagateup=False):
        """Erase all stored data flagged as experimental ('E')"""
        # check if model has been loaded
        if self.model_namespace is None:
            self.load()
        for key in list(self.model_namespace.stored.keys()):
            if self.model_namespace.stored[key][1] == 'E':
                self.model_namespace.stored.pop(key)
        if propagateup:
            for qsar in self.super_models:
                qsar.erase_experimental_stored(propagateup=propagateup)

    def erase_all_stored(self, propagatedown=False, propagateup=False):
        """Erase all stored data, so only predicted values are returned"""
        # check if model has been loaded
        if self.model_namespace is None:
            self.load()
        self.model_namespace.stored = {}
        if propagateup:
            for qsar in self.super_models:
                qsar.erase_all_stored(propagateup=propagateup)

    def apply_model(self, solutes=tuple(), solvents=tuple(), components=tuple(), solutef=tuple(), solventf=tuple(), componentf=tuple()):
        """Take openbabel mol in a list, apply the QSAR and return the results"""
        # check if model has been loaded
        if self.model_namespace is None:
            self.load()
        # first test if this molecule is stored and pass stored results if so
        if solutes[0].normsmiles in self.model_namespace.stored:
            return self.model_namespace.stored[solutes[0].normsmiles]
        # assert that there is the correct number of solutes and solvents
        # current assumption for QSARs is that only one chemical is handled at a time,
        # designated as a solute regardless of property type
        if len(solutes) < self.model_namespace.chemical_inputs['solute min'] or len(solutes) > self.model_namespace.chemical_inputs['solute max'] or \
                len(solvents) < self.model_namespace.chemical_inputs['solvent min'] or len(solvents) > self.model_namespace.chemical_inputs['solvent max'] or \
                len(components) < self.model_namespace.chemical_inputs['component min'] or len(components) > self.model_namespace.chemical_inputs['component max'] or \
                len(solutes)+len(solvents)+len(components) < self.model_namespace.chemical_inputs['total min'] or \
                len(solutes)+len(solvents)+len(components) > self.model_namespace.chemical_inputs['total max']:
            return np.nan, np.nan, np.nan, 'chemical input error: mixture specification not allowed', '', '', ''
        # add or delete hydrogens depending on model
        if self.model_namespace.molecule_format == 'old_format':
            solutes[0].AddHydrogens()
        else:
            solutes[0].DeleteHydrogens()
        # get fragment counts for MLR
        fragment_counts = []
        if self.model_namespace.model_type == 'MLR':
            for smarts in self.model_namespace.smartslist:
                if smarts == 'intercept':
                    fragment_counts.append(1)
                elif smarts == 'sssr':
                    fragment_counts.append(len(solutes[0].GetSSSR()))
                elif smarts == 'MW':
                    fragment_counts.append(solutes[0].GetMolWt())
                else:
                    smarts.Match(solutes[0])
                    fragment_counts.append(len(smarts.GetUMapList()))
        # get fragment counts for MLRX
        elif self.model_namespace.model_type == 'MLRX':
            matchedatoms = set()
            for smarts in self.model_namespace.smartslist:
                if smarts == 'intercept':
                    fragment_counts.append(1)
                elif smarts == 'sssr':
                    fragment_counts.append(len(solutes[0].GetSSSR()))
                elif smarts == 'MW':
                    fragment_counts.append(solutes[0].GetMolWt())
                else:
                    smarts.Match(solutes[0])
                    matchlist = smarts.GetUMapList()
                    matchcount = 0
                    for match in matchlist:
                        if len(set(match).intersection(matchedatoms)) == 0:
                            matchcount += 1
                            matchedatoms.update(set(match))
                    fragment_counts.append(matchcount)
        # get fragment counts for MLRA
        elif self.model_namespace.model_type == 'MLRA':
            for smarts in self.model_namespace.smartslist:
                if smarts == 'intercept':
                    fragment_counts.append(1)
                elif smarts == 'sssr':
                    if len(solutes[0].GetSSSR()) > 0:
                        fragment_counts.append(1)
                    else:
                        fragment_counts.append(0)
                elif smarts == 'MW':
                    fragment_counts.append(solutes[0].GetMolWt())
                else:
                    smarts.Match(solutes[0])
                    if len(smarts.GetUMapList()) > 0:
                        fragment_counts.append(1)
                    else:
                        fragment_counts.append(0)
        fragment_counts = np.array(fragment_counts)
        # apply multiple linear regression qsar
        if self.model_namespace.model_type in ('MLR', 'MLRX', 'MLRA'):
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
                    topgroup.append((fragsim, 1 - self.model_namespace.datalist['value_similarity'][t]))
                    topgroup.sort(reverse=True)
                    if len(topgroup) > topn:
                        for i in reversed(range(topn, len(topgroup))):
                            if topgroup[i][0] != topgroup[topn - 1][0]:
                                topgroup.pop(i)
                    mintop = topgroup[-1][0]
                css = 1
                for i in range(topn):
                    css *= (topgroup[i][0] * (1 - topgroup[i][1])) ** 0.5
                css = css ** (1 / float(topn))
                # calculate leverage
                if self.model_namespace.intercept:
                    leverage = np.matmul(np.matmul(fragment_counts[1:], self.model_namespace.xtxi),
                                         fragment_counts[1:].T)
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
                    pattern1.Match(solutes[0])
                    pattern2.Match(solutes[0])
                    if len(pattern1.GetUMapList()) != len(pattern2.GetUMapList()):
                        violations.append(description)
                if len(violations) > 0:
                    ul = 5
                    error = self.model_namespace.warn_5_error
                    for v in violations:
                        note.append(v)
                # check bounds
                if self.model_namespace.lower_bound and prediction < self.model_namespace.min_train:
                    note.append('prediction ({}) less than smallest value in training set, original aggregate UL: {}'.format(round(prediction, self.model_namespace.round_digits), ul))
                    ul = 6
                    prediction = self.model_namespace.min_train
                if self.model_namespace.upper_bound and prediction > self.model_namespace.max_train:
                    note.append('prediction ({}) greater than largest value in training set, original aggregate UL: {}'.format(round(prediction, self.model_namespace.round_digits), ul))
                    ul = 6
                    prediction = self.model_namespace.max_train
                post_proc_prediction, post_proc_error = self.model_namespace.post_processing(prediction, error)
                self.model_namespace.stored[solutes[0].normsmiles] = (post_proc_prediction, ul, post_proc_error, ', '.join(note), self.model_namespace.citation, self.model_namespace.units, self.model_namespace.endpoint)
                return post_proc_prediction, ul, post_proc_error, ', '.join(note), self.model_namespace.citation, self.model_namespace.units, self.model_namespace.endpoint
            else:
                post_proc_prediction, post_proc_error = self.model_namespace.post_processing(prediction, error)
                self.model_namespace.stored[solutes[0].normsmiles] = (post_proc_prediction, np.nan, post_proc_error, '', self.model_namespace.citation, self.model_namespace.units, self.model_namespace.endpoint)
                return post_proc_prediction, np.nan, post_proc_error, '', self.model_namespace.citation, self.model_namespace.units, self.model_namespace.endpoint


class METAQSARModel:
    """Class that loads a Meta QSAR, which combines data from its dependencies
    to make a new model prediction for IFSMols (subclass of openbabel mols)"""

    def __init__(self, model_module, model_name, version):
        """Save model name and create model namespace"""

        # define Model namespace
        self.model_namespace = None
        self.ismixture = None
        self.model_module = model_module
        self.model_name = model_name
        self.version = version
        self.super_models = []
        self.default_stored = {}

    def __str__(self):
        return self.model_name

    def __repr__(self):
        return self.__str__()

    def load(self):
        """Import the QSAR python module and link the QSAR dependencies"""
        # check if model has been linked
        if self.model_namespace is not None:
            return
        # initiate model namespace
        self.model_namespace = importlib.import_module(self.model_module)
        # check if model is a mixture or not to help format outputs
        if self.model_namespace.chemical_inputs['total min'] <= 1:
            self.ismixture = False
        else:
            self.ismixture = True
        # link models which depend on this one
        # link solute dependencies
        self.model_namespace.solutedependencymodels = {}
        solutedependencies = get_qsar_list(qsarlist=self.model_namespace.solute_dependencies_list)
        for qsar in solutedependencies:
            self.model_namespace.solutedependencymodels[qsar.model_name] = qsar
            self.model_namespace.solutedependencymodels[qsar.model_name].super_models.append(self)
        # link solvent dependencies
        self.model_namespace.solventdependencymodels = {}
        solventdependencies = get_qsar_list(qsarlist=self.model_namespace.solvent_dependencies_list)
        for qsar in solventdependencies:
            self.model_namespace.solventdependencymodels[qsar.model_name] = qsar
            self.model_namespace.solventdependencymodels[qsar.model_name].super_models.append(self)
        # link component dependencies
        self.model_namespace.componentdependencymodels = {}
        componentdependencies = get_qsar_list(qsarlist=self.model_namespace.component_dependencies_list)
        for qsar in componentdependencies:
            self.model_namespace.componentdependencymodels[qsar.model_name] = qsar
            self.model_namespace.componentdependencymodels[qsar.model_name].super_models.append(self)
        # backup the stored data for reset and restore
        self.default_stored = self.model_namespace.stored.copy()

    def set_stored(self, normsmiles, value, ul='U', error=np.nan, ulnote='user value', citation='user value', units=None, endpoint=None):
        """Add a stored value with user-entered data"""
        # check if model has been loaded
        if self.model_namespace is None:
            self.load()
        # do not allow user-entered data for mixtures
        if self.model_namespace.chemical_inputs['total max'] > 1:
            return
        # check that UL is one of the permitted values
        if type(ul) is not int and type(ul) is not str and not np.isnan(ul):
            raise ValueError('User-defined UL must be E, U, int 0-6 or combination thereof')
        elif type(ul) is int and ul not in [0, 1, 2, 3, 4, 5, 6]:
            raise ValueError('User-defined UL must be E, U, int 0-6 or combination thereof')
        elif type(ul) is str:
            counts = []
            for c in ['E', 'U', '0', '1', '2', '3', '4', '5', '6']:
                counts.append(ul.count(c))
                if counts[-1] > 1:
                    raise ValueError('User-defined UL must be E, U, int 0-6 or combination thereof')
            if len(ul) - sum(counts) > 0:
                print(ul, len(ul), counts, sum(counts))
                raise ValueError('User-defined UL must be E, U, int 0-6 or combination thereof')
        # set units and endpoint from model namespace if not provided with input
        if units is None:
            units = self.model_namespace.units
        if endpoint is None:
            endpoint = self.model_namespace.endpoint
        # add user value to stored data
        self.model_namespace.stored[normsmiles] = (value, ul, error, ulnote, citation, units, endpoint)

    def load_stored(self, normsmiles, propagatedown=False, propagateup=False):
        """Remove a specific stored value"""
        # check if model has been loaded
        if self.model_namespace is None:
            self.load()
        if normsmiles in self.default_stored:
            self.model_namespace.stored[normsmiles] = self.default_stored[normsmiles]
        if propagatedown:
            for qsar in self.model_namespace.solutedependencymodels.values():
                qsar.load_stored(normsmiles, propagatedown=propagatedown)
            for qsar in self.model_namespace.componentdependencymodels.values():
                qsar.load_stored(normsmiles, propagatedown=propagatedown)
            for qsar in self.model_namespace.solventdependencymodels.values():
                qsar.load_stored(normsmiles, propagatedown=propagatedown)
        if propagateup:
            for qsar in self.super_models:
                qsar.load_stored(normsmiles, propagateup=propagateup)

    def remove_stored(self, normsmiles, propagatedown=False, propagateup=False):
        """Remove a specific stored value"""
        # check if model has been loaded
        if self.model_namespace is None:
            self.load()
        if normsmiles in self.model_namespace.stored:
            self.model_namespace.stored.pop(normsmiles)
        if propagatedown:
            for qsar in self.model_namespace.solutedependencymodels.values():
                qsar.remove_stored(normsmiles, propagatedown=propagatedown)
            for qsar in self.model_namespace.componentdependencymodels.values():
                qsar.remove_stored(normsmiles, propagatedown=propagatedown)
            for qsar in self.model_namespace.solventdependencymodels.values():
                qsar.remove_stored(normsmiles, propagatedown=propagatedown)
        if propagateup:
            for qsar in self.super_models:
                qsar.remove_stored(normsmiles, propagateup=propagateup)

    def reset_stored(self, propagatedown=False, propagateup=False):
        """Reset the stored data to default, erasing all predicted and user data"""
        # check if model has been loaded
        if self.model_namespace is None:
            self.load()
        self.model_namespace.stored = self.default_stored.copy()
        if propagatedown:
            for qsar in self.model_namespace.solutedependencymodels.values():
                qsar.reset_stored(propagatedown=propagatedown)
            for qsar in self.model_namespace.componentdependencymodels.values():
                qsar.reset_stored(propagatedown=propagatedown)
            for qsar in self.model_namespace.solventdependencymodels.values():
                qsar.reset_stored(propagatedown=propagatedown)
        if propagateup:
            for qsar in self.super_models:
                qsar.reset_stored(propagateup=propagateup)

    def restore_stored(self, propagatedown=False, propagateup=False):
        """Reset the stored data to default, preserving predicted and user data not in the default data"""
        # check if model has been loaded
        if self.model_namespace is None:
            self.load()
        self.model_namespace.stored.update(self.default_stored)
        if propagatedown:
            for qsar in self.model_namespace.solutedependencymodels.values():
                qsar.restore_stored(propagatedown=propagatedown)
            for qsar in self.model_namespace.componentdependencymodels.values():
                qsar.restore_stored(propagatedown=propagatedown)
            for qsar in self.model_namespace.solventdependencymodels.values():
                qsar.restore_stored(propagatedown=propagatedown)
        if propagateup:
            for qsar in self.super_models:
                qsar.restore_stored(propagateup=propagateup)

    def erase_user_stored(self, propagatedown=False, propagateup=False):
        """Erase all stored data flagged as user ('U')"""
        # check if model has been loaded
        if self.model_namespace is None:
            self.load()
        for key in list(self.model_namespace.stored.keys()):
            if self.model_namespace.stored[key][1] == 'U':
                self.model_namespace.stored.pop(key)
        if propagatedown:
            for qsar in self.model_namespace.solutedependencymodels.values():
                qsar.erase_user_stored(propagatedown=propagatedown)
            for qsar in self.model_namespace.componentdependencymodels.values():
                qsar.erase_user_stored(propagatedown=propagatedown)
            for qsar in self.model_namespace.solventdependencymodels.values():
                qsar.erase_user_stored(propagatedown=propagatedown)
        if propagateup:
            for qsar in self.super_models:
                qsar.erase_user_stored(propagateup=propagateup)

    def erase_experimental_stored(self, propagatedown=False, propagateup=False):
        """Erase all stored data flagged as experimental ('E')"""
        # check if model has been loaded
        if self.model_namespace is None:
            self.load()
        for key in list(self.model_namespace.stored.keys()):
            if self.model_namespace.stored[key][1] == 'E':
                self.model_namespace.stored.pop(key)
        if propagatedown:
            for qsar in self.model_namespace.solutedependencymodels.values():
                qsar.erase_experimental_stored(propagatedown=propagatedown)
            for qsar in self.model_namespace.componentdependencymodels.values():
                qsar.erase_experimental_stored(propagatedown=propagatedown)
            for qsar in self.model_namespace.solventdependencymodels.values():
                qsar.erase_experimental_stored(propagatedown=propagatedown)
        if propagateup:
            for qsar in self.super_models:
                qsar.erase_experimental_stored(propagateup=propagateup)

    def erase_all_stored(self, propagatedown=False, propagateup=False):
        """Erase all stored data, so only predicted values are returned"""
        # check if model has been loaded
        if self.model_namespace is None:
            self.load()
        self.model_namespace.stored = {}
        if propagatedown:
            for qsar in self.model_namespace.solutedependencymodels.values():
                qsar.erase_all_stored(propagatedown=propagatedown)
            for qsar in self.model_namespace.componentdependencymodels.values():
                qsar.erase_all_stored(propagatedown=propagatedown)
            for qsar in self.model_namespace.solventdependencymodels.values():
                qsar.erase_all_stored(propagatedown=propagatedown)
        if propagateup:
            for qsar in self.super_models:
                qsar.erase_all_stored(propagateup=propagateup)

    def apply_model(self, solutes=tuple(), solvents=tuple(), components=tuple(), solutef=tuple(), solventf=tuple(), componentf=tuple()):
        """Take openbabel mol(s) in lists of solutes and solvents, apply the Meta QSAR and return the results"""
        # check if model has been loaded and dependencies linked
        if self.model_namespace is None:
            self.load()
        # assert that there is the correct number of solutes and solvents
        if len(solutes) < self.model_namespace.chemical_inputs['solute min'] or len(solutes) > self.model_namespace.chemical_inputs['solute max'] or \
                len(solvents) < self.model_namespace.chemical_inputs['solvent min'] or len(solvents) > self.model_namespace.chemical_inputs['solvent max'] or \
                len(components) < self.model_namespace.chemical_inputs['component min'] or len(components) > self.model_namespace.chemical_inputs['component max'] or \
                len(solutes)+len(solvents)+len(components) < self.model_namespace.chemical_inputs['total min'] or \
                len(solutes)+len(solvents)+len(components) > self.model_namespace.chemical_inputs['total max']:
            if self.model_namespace.chemical_inputs['total max'] == 1:
                return np.nan, np.nan, np.nan, 'chemical input error: mixture specification not allowed', '', '', ''
            elif self.model_namespace.chemical_inputs['total max'] == 2 and \
                    self.model_namespace.chemical_inputs['solute min'] == 1 and \
                    self.model_namespace.chemical_inputs['solvent min'] == 1:
                return np.nan, np.nan, np.nan, 'chemical input error: mixture specification with one solute and one solvent required', '', '', ''
            elif self.model_namespace.chemical_inputs['total min'] >= 2 and \
                     self.model_namespace.chemical_inputs['component min'] + self.model_namespace.chemical_inputs['solvent min'] >= 2:
                return np.nan, np.nan, np.nan, 'chemical input error: mixture specification with at least two components/solvents required', '', '', ''
            else:
                return np.nan, np.nan, np.nan, 'chemical input error: mixture specification required', '', '', ''
        # check if a value is stored for this solute
        if self.model_namespace.chemical_inputs['total max'] == 1 and solutes[0].normsmiles in self.model_namespace.stored:
            return self.model_namespace.stored[solutes[0].normsmiles]
        # # check if a value is stored for this solute/solvent pair
        # elif self.model_namespace.chemical_inputs['total max'] == 2 and \
        #         self.model_namespace.chemical_inputs['solute min'] == 1 and \
        #         self.model_namespace.chemical_inputs['solvent min'] == 1:
        #     key = '|'.join([solutes[0].normsmiles, solvents[0].normsmiles])
        #     if key in self.model_namespace.stored:
        #         return self.model_namespace.stored[key]
        # # check if this exact mixture has a stored value
        # elif self.model_namespace.chemical_inputs['total min'] >= 2:
        #     keyparts = []
        #     for i in range(len(solutes)):
        #         keyparts.append('solute')
        #         keyparts.append(solutes[i].normsmiles)
        #         keyparts.append(solutef[i][0])
        #         keyparts.append(solutef[i][1])
        #     for i in range(len(solvents)):
        #         keyparts.append('solvent')
        #         keyparts.append(solvents[i].normsmiles)
        #         keyparts.append(solventf[i][0])
        #         keyparts.append(solventf[i][1])
        #     for i in range(len(components)):
        #         keyparts.append('component')
        #         keyparts.append(components[i].normsmiles)
        #         keyparts.append(componentf[i][0])
        #         keyparts.append(componentf[i][1])
        #     key = '|'.join(keyparts)
        #     if key in self.model_namespace.stored:
        #         return self.model_namespace.stored[key]
        # pass solute to solute dependency qsar models to calculate values
        solutedependencies = []
        for s in range(len(solutes)):
            solutedependencies.append({})
            for d, m in self.model_namespace.solutedependencymodels.items():
                solutedependencies[-1][d] = tuple(m.apply_model(solutes=(solutes[s],), solutef=(solutef[s],)))
        # pass solvents to solvent dependency qsar models to calculate values
        solventdependencies = []
        for s in range(len(solvents)):
            solventdependencies.append({})
            for d, m in self.model_namespace.solventdependencymodels.items():
                solventdependencies[-1][d] = tuple(m.apply_model(solutes=(solvents[s],), solutef=(solventf[s],)))
        # pass components to component dependency qsar models to calculate values
        componentdependencies = []
        for c in range(len(components)):
            componentdependencies.append({})
            for d, m in self.model_namespace.componentdependencymodels.items():
                componentdependencies[-1][d] = tuple(m.apply_model(solutes=(components[c],), solutef=(componentf[c],)))
        # generate propagated domain notes for solutes
        propagated_domain_notes = []
        for s in range(len(solutes)):
            domainnotes = []
            for k, v in solutedependencies[s].items():
                if not type(solutedependencies[s][k][1]) is str and np.isnan(solutedependencies[s][k][1]):
                    continue
                domainnotes.append(''.join([k, '=', str(solutedependencies[s][k][1])]))
            if len(domainnotes):
                propagated_domain_notes.append(''.join(['solute {} dependency ULs: '.format(s+1), ', '.join(domainnotes)]))
        # generate propagated domain notes for components
        for c in range(len(componentdependencies)):
            domainnotes = []
            for k, v in componentdependencies[c].items():
                if not type(componentdependencies[c][k][1]) is str and np.isnan(componentdependencies[c][k][1]):
                    continue
                domainnotes.append(''.join([k, '=', str(componentdependencies[c][k][1])]))
            if len(domainnotes):
                propagated_domain_notes.append(''.join(['component {} dependency ULs: '.format(c+1), ', '.join(domainnotes)]))
        # generate propagated domain notes for solvents
        for s in range(len(solvents)):
            domainnotes = []
            for k, v in solventdependencies[s].items():
                if not type(solventdependencies[s][k][1]) is str and np.isnan(solventdependencies[s][k][1]):
                    continue
                domainnotes.append(''.join([k, '=', str(solventdependencies[s][k][1])]))
            if len(domainnotes):
                propagated_domain_notes.append(''.join(['solvent {} dependency ULs: '.format(s+1), ', '.join(domainnotes)]))
        self.model_namespace.propagated_domain_notes = '; '.join(propagated_domain_notes)
        # call the metamodel with the dependency outputs to calculate meta result
        prediction, UL, error, ULnote, citation, units, endpoint = self.model_namespace.calculate(solutedependencies, solventdependencies, componentdependencies, solutef, solventf, componentf)
        # store result for this solute
        if self.model_namespace.chemical_inputs['total max'] == 1:
            self.set_stored(solutes[0].normsmiles, prediction, UL, error, ULnote, citation, units, endpoint)
        # # store result for this solute/solvent pair
        # elif self.model_namespace.chemical_inputs['total max'] == 2 and \
        #         self.model_namespace.chemical_inputs['solute min'] == 1 and \
        #         self.model_namespace.chemical_inputs['solvent min'] == 1:
        #     key = '|'.join([solutes[0].normsmiles, solvents[0].normsmiles])
        #     self.model_namespace.stored[key] = (prediction, UL, error, ULnote, citation, units, endpoint)
        # # store result for this exact mixture
        # elif self.model_namespace.chemical_inputs['total min'] >= 2:
        #     keyparts = []
        #     for i in range(len(solutes)):
        #         keyparts.append('solute')
        #         keyparts.append(solutes[i].normsmiles)
        #         keyparts.append(solutef[i][0])
        #         keyparts.append(solutef[i][1])
        #     for i in range(len(solvents)):
        #         keyparts.append('solvent')
        #         keyparts.append(solvents[i].normsmiles)
        #         keyparts.append(solventf[i][0])
        #         keyparts.append(solventf[i][1])
        #     for i in range(len(components)):
        #         keyparts.append('component')
        #         keyparts.append(components[i].normsmiles)
        #         keyparts.append(componentf[i][0])
        #         keyparts.append(componentf[i][1])
        #     key = '|'.join(keyparts)
        #     self.model_namespace.stored[key] = (prediction, UL, error, ULnote, citation, units, endpoint)
        # return result
        return prediction, UL, error, ULnote, citation, units, endpoint


# instantiate qsar models
fhlb = QSARModel('ifsqsar.models.ifs_qsar_fhlb_linr', 'fhlb', 1)
hhlb = QSARModel('ifsqsar.models.ifs_qsar_hhlb_linr', 'hhlb', 1)
hhlt = QSARModel('ifsqsar.models.ifs_qsar_hhlt_linr', 'hhlt', 1)
biowin3usmmlrx = QSARModel('ifsqsar.models.other_qsar_biowin3_usm_mlrx', 'biowin3usmmlrx', 1)
biowin3usmmlra = QSARModel('ifsqsar.models.other_qsar_biowin3_usm_mlra', 'biowin3usmmlra', 1)
biowin4psmmlrx = QSARModel('ifsqsar.models.other_qsar_biowin4_psm_mlrx', 'biowin4psmmlrx', 1)
biowin4psmmlra = QSARModel('ifsqsar.models.other_qsar_biowin4_psm_mlra', 'biowin4psmmlra', 1)
HLbiodeg = METAQSARModel('ifsqsar.models.meta_qsar_hlbiodeg', 'HLbiodeg', 1)
dsm = QSARModel('ifsqsar.models.ifs_qsar_dsm_linr', 'dsm', 1)
tm = QSARModel('ifsqsar.models.ifs_qsar_tm_linr', 'tm', 1)
tmpplfer = METAQSARModel('ifsqsar.models.meta_qsar_tm_pplfer', 'tmpplfer', 1)
tmconsensus = METAQSARModel('ifsqsar.models.meta_qsar_tm_consensus', 'tmconsensus', 1)
tbpplfer = METAQSARModel('ifsqsar.models.meta_qsar_tb_pplfer', 'tbpplfer', 1)
Ev1 = QSARModel('ifsqsar.models.ifs_qsar_ADB_UFZ__E_linr', 'E', 1)
Sv1 = QSARModel('ifsqsar.models.ifs_qsar_ADB_UFZ__S_linr', 'S', 1)
Av1 = QSARModel('ifsqsar.models.ifs_qsar_ADB_UFZ__A_linr', 'A', 1)
Bv1 = QSARModel('ifsqsar.models.ifs_qsar_ADB_UFZ__B_linr', 'B', 1)
Lv1 = QSARModel('ifsqsar.models.ifs_qsar_ADB_UFZ__L_linr', 'L', 1)
Vtd = QSARModel('ifsqsar.models.other_qsar_V', 'V', 1)
Ev2 = QSARModel('ifsqsar.models.ifs_qsar_pplfer_solutes_E_linr', 'E', 2)
Sv2 = QSARModel('ifsqsar.models.ifs_qsar_pplfer_solutes_S_linr', 'S', 2)
Av2 = QSARModel('ifsqsar.models.ifs_qsar_pplfer_solutes_A_linr', 'A', 2)
Bv2 = QSARModel('ifsqsar.models.ifs_qsar_pplfer_solutes_B_linr', 'B', 2)
Lv2 = QSARModel('ifsqsar.models.ifs_qsar_pplfer_solutes_L_linr', 'L', 2)
ssp = QSARModel('ifsqsar.models.ifs_qsar_pplfer_system_2_s_linr', 's', 1)
asp = QSARModel('ifsqsar.models.ifs_qsar_pplfer_system_2_a_linr', 'a', 1)
bsp = QSARModel('ifsqsar.models.ifs_qsar_pplfer_system_2_b_linr', 'b', 1)
vsp = QSARModel('ifsqsar.models.ifs_qsar_pplfer_system_2_v_linr', 'v', 1)
lsp = QSARModel('ifsqsar.models.ifs_qsar_pplfer_system_2_l_linr', 'l', 1)
csp = QSARModel('ifsqsar.models.ifs_qsar_pplfer_system_2_c_linr', 'c', 1)
logKow = METAQSARModel('ifsqsar.models.meta_qsar_logkow_pplfer', 'logKow', 1)
logKowdry = METAQSARModel('ifsqsar.models.meta_qsar_logkowdry_pplfer', 'logKowdry', 1)
logKoa = METAQSARModel('ifsqsar.models.meta_qsar_logkoa_pplfer', 'logKoa', 1)
logKaw = METAQSARModel('ifsqsar.models.meta_qsar_logkaw_pplfer', 'logKaw', 1)
logKoo = METAQSARModel('ifsqsar.models.meta_qsar_logkowod_pplfer', 'logKoo', 1)
logVPliquid = METAQSARModel('ifsqsar.models.meta_qsar_logVP_liquid', 'logVPliquid', 1)
logSwliquid = METAQSARModel('ifsqsar.models.meta_qsar_logSw_liquid', 'logSwliquid', 1)
logSoliquid = METAQSARModel('ifsqsar.models.meta_qsar_logSo_liquid', 'logSoliquid', 1)
logSowetliquid = METAQSARModel('ifsqsar.models.meta_qsar_logSowet_liquid', 'logSowetliquid', 1)
logVPliquidpplfer = METAQSARModel('ifsqsar.models.meta_qsar_logVP_liquid_pplfer', 'logVPliquid', 2)
logSwliquidpplfer = METAQSARModel('ifsqsar.models.meta_qsar_logSw_liquid_pplfer', 'logSwliquid', 2)
logKsa = METAQSARModel('ifsqsar.models.meta_qsar_logksa_pplfer', 'logKsa', 1)
MVmlrx = QSARModel('ifsqsar.models.other_qsar_MV_mlrx', 'MVmlrx', 1)
MVmlr = QSARModel('ifsqsar.models.other_qsar_MV_mlr', 'MVmlr', 1)
MVmlrRings = QSARModel('ifsqsar.models.other_qsar_MV_mlrRings', 'MVmlrRings', 1)
MVliqcorr = QSARModel('ifsqsar.models.other_qsar_MV_liq_corr', 'MVliqcorr', 1)
MVsolid = METAQSARModel('ifsqsar.models.meta_qsar_MV_solid', 'MVsolid', 1)
MVliquid = METAQSARModel('ifsqsar.models.meta_qsar_MV_liquid', 'MVliquid', 1)
MVsolidV = METAQSARModel('ifsqsar.models.meta_qsar_MV_solid_V', 'MVsolid', 1)
MVliquidV = QSARModel('ifsqsar.models.other_qsar_MV_liquid_V', 'MVliquid', 1)
densitysolid = METAQSARModel('ifsqsar.models.meta_qsar_density_solid', 'densitysolid', 1)
densityliquid = METAQSARModel('ifsqsar.models.meta_qsar_density_liquid', 'densityliquid', 1)
MW = QSARModel('ifsqsar.models.other_qsar_MW', 'MW', 1)
state = METAQSARModel('ifsqsar.models.meta_qsar_state', 'state', 1)


def get_qsar_list(qsarlist=None, versionlist=None):
    """Main interface for getting lists of QSARs which can be filtered by optional selection criteria.

    Arguments:
        qsarlist -- a list of models names (strings) from below (default is all current models)
            "fhlb" -- biotransformation half-life in fish
            "hhlb" -- biotransformation half-life in humans
            "hhlt" -- total elimination half-life in humans
            "HLbiodeg" -- aqueous aerobic biodegradation half-life
            "dsm" -- entropy of melting (aka entropy of fusion)
            "tm" -- melting point (QSPR)
            "tmpplfer" -- melting point (PPLFER)
            "tmconsensus" -- melting point (mean of QSPR+PPLFER)
            "E","S","A","B","L","V" -- Abraham PPLFER solute descriptors
            "s","a","b","v","l","c" -- Abraham/Goss PPLFER system parameters
            "logKow","logKowdry","logKoa","logKaw" -- partition coeff. of octanol, water, and air
            "logKoo" -- octanol wet-dry conversion factor
            "logVPliquid" -- Vapor Pressure of liquids (PPLFER)
            "logSwliquid" -- Water Solubilility of liquids (PPLFER)
            "logSoliquid" -- Octanol Solubility of liquids (PPLFER)
            "MVliquid, MVsolid" -- Molar Volume of liquids and solids
            "densityliquid, densitysolid" -- density of liquids and solids
            "MW" -- Molecular Weight
            "state" -- chemical state (gas, liquid or solid)
            "logKsa" -- solvent-air partitioning (solute-solvent pair, needs mixture spec. SMILES)
        versionlist -- list of version numbers of QSARs as integers, by default the most recent version is returned
    """

    # restrictions based on what is implemented
    # when requesting old model versions a list of qsars must also be passed
    if versionlist is not None and qsarlist is None:
        raise RuntimeError('When requesting old QSAR versions from get_qsar_list a list of'
                           ' QSARs and a list of versions both must be passed.')
    # when requesting old model versions the list of qsars and versions must be the same shape
    elif versionlist is not None and qsarlist is not None and len(versionlist) != len(qsarlist):
        raise RuntimeError('When calling get_qsar_list qsarlist and versionlist must the same length.')
    # when requesting old versions you cannot request both the old and new versions at the same time
    elif versionlist is not None and len(qsarlist) != len(set(qsarlist)):
        raise RuntimeError('When calling get_qsar_list you cannot request multiple versions of the same QSAR.')

    # if no qsar list is passed set the full default list
    if qsarlist is None:
        qsarlist = ['fhlb', 'hhlb', 'hhlt', 'HLbiodeg',
                    'dsm', 'tmconsensus', 'tbpplfer',
                    'logKow', 'logKoa', 'logKaw', 'logKoo',
                    'logVPliquid', 'logSwliquid', 'logSoliquid',
                    'MVliquid', 'densityliquid', 'MW',
                    'state',
                    'E', 'S', 'A', 'B', 'V', 'L',
                    's', 'a', 'b', 'v', 'l', 'c',
                    'logKsa']
    returnlist = []
    # decide if old versions are included in parse list
    currentqsarversions = [fhlb, hhlb, hhlt, biowin3usmmlrx, biowin3usmmlra, biowin4psmmlrx, biowin4psmmlra, HLbiodeg,
                           dsm, tm, tmpplfer, tbpplfer, tmconsensus, state,
                           Ev2, Sv2, Av2, Bv2, Lv2, Vtd, ssp, asp, bsp, vsp, lsp, csp,
                           logKow, logKowdry, logKoa, logKaw, logKoo,
                           logVPliquidpplfer, logSwliquidpplfer, logSoliquid, logSowetliquid,
                           MVsolidV, MVliquidV, densitysolid, densityliquid, MW,
                           logKsa]
    if versionlist is None:
        oldqsarversions = []
    else:
        oldqsarversions = [Ev1, Sv1, Av1, Bv1, Lv1,
                           logVPliquid, logSwliquid,
                           MVmlrx, MVmlr, MVmlrRings, MVliqcorr, MVsolid, MVliquid]
    # parse through all the loaded QSARs and check if they should be added to returnlist
    for i in currentqsarversions + oldqsarversions:
        # first check if QSAR is in passed list of names
        if qsarlist is None or (qsarlist is not None and i.model_name in qsarlist):
            # check if QSAR has version in passed specification, take most recent version by default
            if versionlist is None:
                returnlist.append(i)
            else:
                version = versionlist[qsarlist.index(i.model_name)]
                if i.version == version:
                    returnlist.append(i)

    # sort the list to be in the same order as the input list
    returnlist.sort(key=lambda x: qsarlist.index(x.model_name))

    return returnlist

