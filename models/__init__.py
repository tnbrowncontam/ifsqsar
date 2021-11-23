from openbabel import openbabel as ob
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
    return 1. / (np.exp(-x * 358. / 23. + 111 * np.arctan(x * 37. / 294.)) + 1.)


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
    """Class that loads an arbitrary Model from a Model file and applies it
    to any molecules passed to it as an openbabel mol"""

    def __init__(self, model_module, model_name):
        """Save model name and create model namespace."""
        # define Model namespace
        self.model_namespace = None
        self.model_module = model_module
        self.model_name = model_name
        self.last_molecule = None
        self.last_result = (None, None, None, None)

    def load(self):
        """Import the model module."""
        if self.model_namespace is not None:
            return
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

    def apply_model(self, solutes, solvents=tuple()):
        """Take an openbabel molecule, apply the QSARModel and return the result."""
        # first test if this molecule is last molecule and pass last results if so
        if solutes is self.last_molecule:
            return self.last_result[0], self.last_result[1], self.last_result[2], self.last_result[3]
        # check if model has been loaded
        if self.model_namespace is None:
            self.load()
        # only one chemical handled at a time, designated as a solute regardless of property type
        assert len(solutes) == self.model_namespace.components['solute'] and len(solvents) == self.model_namespace.components['solvent']
        # add or delete hydrogens depending on model
        if self.model_namespace.molecule_format == 'old_format':
            solutes[0].AddHydrogens()
        else:
            solutes[0].DeleteHydrogens()
        # get fragment counts
        fragment_counts = []
        for smarts in self.model_namespace.smartslist:
            if smarts == 'intercept':
                fragment_counts.append(1)
            else:
                smarts.Match(solutes[0])
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
                    ul = 6
                    note.append('prediction (' + str(round(prediction,
                                                           self.model_namespace.round_digits)) + ') less than smallest value in training set')
                    prediction = self.model_namespace.min_train
                if self.model_namespace.upper_bound and prediction > self.model_namespace.max_train:
                    ul = 6
                    note.append('prediction (' + str(round(prediction,
                                                           self.model_namespace.round_digits)) + ') greater than largest value in training set')
                    prediction = self.model_namespace.max_train
                post_proc_prediction, post_proc_error = self.model_namespace.post_processing(prediction, error)
                self.last_molecule = solutes
                self.last_result = (post_proc_prediction, ul, post_proc_error, ', '.join(note))
                return post_proc_prediction, ul, post_proc_error, ', '.join(note)
            else:
                post_proc_prediction, post_proc_error = self.model_namespace.post_processing(prediction, error)
                self.last_molecule = solutes
                self.last_result = (post_proc_prediction, np.nan, post_proc_error, '')
                return post_proc_prediction, np.nan, post_proc_error, ''


class METAQSARModel:
    """Class that loads an arbitrary Metamodel from a Metamodel file and applies it
    to any molecules passed to it as an openbabel mol"""

    def __init__(self, model_module, model_name):
        """Save model name and create model namespace."""

        # define Model namespace
        self.model_namespace = None
        self.model_module = model_module
        self.model_name = model_name
        self.last_solute = None
        self.last_solvent = None
        self.last_result = (None, None, None, None)

    def load(self):
        """Import the model module and link the QSAR dependencies."""
        # check if model has been linked
        if self.model_namespace is not None:
            return
        # initiate model namespace
        self.model_namespace = importlib.import_module(self.model_module)
        # link solute dependencies
        self.model_namespace.solutedependencymodels = {}
        solutedependencies = get_qsar_list(qsarlist=self.model_namespace.solute_dependencies_list)
        for qsar in solutedependencies:
            self.model_namespace.solutedependencymodels[qsar.model_name] = qsar
        # link solvent dependencies
        self.model_namespace.solventdependencymodels = {}
        solventdependencies = get_qsar_list(qsarlist=self.model_namespace.solvent_dependencies_list)
        for qsar in solventdependencies:
            self.model_namespace.solventdependencymodels[qsar.model_name] = qsar

    def apply_model(self, solutes=tuple(), solvents=tuple()):
        """Take an openbabel molecule, apply the METAQSARModel and return the result."""
        # first test if this solute and solvent are the same as last calculation and pass last results if so
        if solutes is self.last_solute and solvents is self.last_solvent:
            return self.last_result[0], self.last_result[1], self.last_result[2], self.last_result[3]
        # check if model has been loaded and dependencies linked
        if self.model_namespace is None:
            self.load()
        # assert that there is the correct number of solutes and solvents
        assert len(solutes) == self.model_namespace.components['solute'] and len(solvents) == self.model_namespace.components['solvent']
        # pass solute to solute dependency qsar models to calculate values
        solutedependencies = {}
        for d, m in self.model_namespace.solutedependencymodels.items():
            solutedependencies[d] = tuple(m.apply_model(solutes))
        # pass solvents to solvent dependency qsar models to calculate values
        solventdependencies = {}
        for d, m in self.model_namespace.solventdependencymodels.items():
            solutedependencies[d] = tuple(m.apply_model(solvents))
        # call the metamodel with the dependency outputs to calculate meta result
        prediction, UL, error, ULnote = self.model_namespace.calculate(solutedependencies, solventdependencies)
        # return result
        self.last_solute = solutes
        self.last_solvent = solvents
        self.last_result = (prediction, UL, error, ULnote)
        return prediction, UL, error, ULnote


# instantiate qsar models
fhlb = QSARModel('ifsqsar.models.ifs_qsar_fhlb_linr', 'fhlb')
hhlb = QSARModel('ifsqsar.models.ifs_qsar_hhlb_linr', 'hhlb')
hhlt = QSARModel('ifsqsar.models.ifs_qsar_hhlt_linr', 'hhlt')
dsm = QSARModel('ifsqsar.models.ifs_qsar_dsm_linr', 'dsm')
tm = QSARModel('ifsqsar.models.ifs_qsar_tm_linr', 'tm')
Ev1 = QSARModel('ifsqsar.models.ifs_qsar_ADB_UFZ__E_linr', 'E')
Sv1 = QSARModel('ifsqsar.models.ifs_qsar_ADB_UFZ__S_linr', 'S')
Av1 = QSARModel('ifsqsar.models.ifs_qsar_ADB_UFZ__A_linr', 'A')
Bv1 = QSARModel('ifsqsar.models.ifs_qsar_ADB_UFZ__B_linr', 'B')
Lv1 = QSARModel('ifsqsar.models.ifs_qsar_ADB_UFZ__L_linr', 'L')
Vtd = QSARModel('ifsqsar.models.ifs_qsar_V', 'V')
Ev2 = QSARModel('ifsqsar.models.ifs_qsar_pplfer_solutes_E_linr', 'E')
Sv2 = QSARModel('ifsqsar.models.ifs_qsar_pplfer_solutes_S_linr', 'S')
Av2 = QSARModel('ifsqsar.models.ifs_qsar_pplfer_solutes_A_linr', 'A')
Bv2 = QSARModel('ifsqsar.models.ifs_qsar_pplfer_solutes_B_linr', 'B')
Lv2 = QSARModel('ifsqsar.models.ifs_qsar_pplfer_solutes_L_linr', 'L')
svd = QSARModel('ifsqsar.models.ifs_qsar_pplfer_system_2_s_linr', 's')
avd = QSARModel('ifsqsar.models.ifs_qsar_pplfer_system_2_s_linr', 'a')
bvd = QSARModel('ifsqsar.models.ifs_qsar_pplfer_system_2_s_linr', 'b')
vvd = QSARModel('ifsqsar.models.ifs_qsar_pplfer_system_2_s_linr', 'v')
lvd = QSARModel('ifsqsar.models.ifs_qsar_pplfer_system_2_s_linr', 'l')
cvd = QSARModel('ifsqsar.models.ifs_qsar_pplfer_system_2_s_linr', 'c')
logKow = METAQSARModel('ifsqsar.models.meta_qsar_logkow_pplfer', 'logKow')


def get_qsar_list(qsarlist=None, components=None, version=None):
    """function for getting lists of QSARs meeting selection criteria"""
    returnlist = []
    # calculate number of solutes and solvents in passed components specification
    solutenumber = 0
    solventnumber = 0
    if components is not None:
        if 'solute' in components:
            solutenumber = components['solute']
        if 'component' in components:
            solutenumber = components['solute']
            solventnumber = components['solvent']
        if 'solvent' in components:
            solventnumber = components['solvent']
    # decide if old versions are included in parse list
    currentqsarversions = [fhlb, hhlb, hhlt, dsm, tm, Ev2, Sv2, Av2, Bv2, Lv2, Vtd, logKow]
    if version is None:
        oldqsarversions = []
    else:
        oldqsarversions = [Ev1, Sv1, Av1, Bv1, Lv1]
    # parse through all the loaded QSARs and check if they should be added to returnlist
    for i in currentqsarversions + oldqsarversions:
        # first check if QSAR is in passed list of names
        if qsarlist is None or (qsarlist is not None and i.model_name in qsarlist):
            # check if the QSAR handles the number of solutes and solvents passed
            if components is None or \
                    ((i.model_namespace.components['solute'] is None or solutenumber == i.model_namespace.components['solute']) and
                     (i.model_namespace.components['solvent'] is None or solventnumber == i.model_namespace.components['solvent'])):
                # check if QSAR has version in passed specification, take most recent version by default
                if version is None or (version is not None and i.model_namespace.version[0] == version):
                    returnlist.append(i)
    return returnlist

