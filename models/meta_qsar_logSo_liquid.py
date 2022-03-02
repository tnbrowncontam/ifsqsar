"""Meta QSAR for logSoliquid (dry solvent)"""
import numpy as np
value_names = ('logSoliquid',)
version = 1
endpoint = 'Log of solubility in octanol for liquid or super-cooled liquid solute'
citation = 'Solute descriptors and VP: '\
           'Brown, T. N.; '\
           'Development of Iterative Fragment Selection (IFS) QSPRs for Poly-Parameter Linear Free Energy '\
           'Relationship (PPLFER) Solute Descriptors and System Parameters. '\
           'J Solution Chem 2021, In Review. '\
           'PPLFER Equation for thermodynamic cycle: '\
           'Brown, T. N., '\
           'Empirical Regressions between System Parameters and Solute Descriptors of Polyparameter Linear Free '\
           'Energy Relationships (PPLFERs) for Predicting Solvent-Air Partitioning. '\
           'Fluid Phase Equilibria 2021, 113035.'
round_digits = 2
units = 'log mg/L'
components = {'solute': 1, 'solvent': 0}
solute_dependencies_list = ['logVPliquid', 'logKoa', 'MW']
solvent_dependencies_list = []
propagated_domain_notes = ''
smiles_flag = 'neutrals'


def calculate(solutedependencies, solventdependencies):
    # convert VP to Octanol Solubility by thermodynamic cycle
    logSo = solutedependencies['logVPliquid'][0] - np.log10(8.31446261815324) - np.log10(293.15) + solutedependencies['logKoa'][0]
    # convert from mol/m3 to mg/L
    logSo = logSo + np.log10((1000 * solutedependencies['MW'][0]) * (1/1000))
    # get average UL
    logSoUL = int(np.ceil(((solutedependencies['logVPliquid'][1]**2 + solutedependencies['logKoa'][1]**2)/2)**0.5))
    # calculate total error
    logSoerr = (solutedependencies['logVPliquid'][2]**2 + solutedependencies['logKoa'][2]**2)**0.5
    # set domain notes
    domainnotes = [propagated_domain_notes]
    if logSoUL <= 1:
        domainnotes.append('aggregate solute descriptor UL is in the AD')
    else:
        domainnotes.append('aggregate solute descriptor UL is out of the AD')

    return round(logSo, round_digits), logSoUL, round(logSoerr, round_digits), '; '.join(domainnotes), citation, units, endpoint

