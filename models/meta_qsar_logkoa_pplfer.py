"""Meta QSAR for logKoa (dry solvent)"""
import numpy as np
value_names = ('logKoa',)
version = 1
citation = 'Solute descriptors: '\
           'Brown, T. N.; '\
           'Development of Iterative Fragment Selection (IFS) QSPRs for Poly-Parameter Linear Free Energy '\
           'Relationship (PPLFER) Solute Descriptors and System Parameters. '\
           'J Solution Chem 2021, In Review. '\
           'PPLFER Equation: '\
           'Brown, T. N., '\
           'Empirical Regressions between System Parameters and Solute Descriptors of Polyparameter Linear Free '\
           'Energy Relationships (PPLFERs) for Predicting Solvent-Air Partitioning. '\
           'Fluid Phase Equilibria 2021, 113035.'
round_digits = 2
units = 'log L[a]/L[o]'
components = {'solute': 1, 'solvent': 0}
solute_dependencies_list = ['S', 'A', 'B', 'V', 'L']
solvent_dependencies_list = []
propagated_domain_notes = ''
smiles_flag = 'neutrals'


def calculate(solutedependencies, solventdependencies):
    # calculate log Koa with pplfer equation from Brown 2021
    logKoa = 0.69 * solutedependencies['S'][0] \
             + 3.56 * solutedependencies['A'][0] \
             + 0.73 * solutedependencies['B'][0] \
             + 0.52 * solutedependencies['V'][0] \
             + 0.79 * solutedependencies['L'][0] \
             -0.26
    logKoaUL = 0
    for sltdes in ['S', 'A', 'B', 'L']:
        if solutedependencies[sltdes][1] < 4:
            logKoaUL += solutedependencies[sltdes][1]**2
        elif solutedependencies[sltdes][1] == 4 and sltdes != 'L':
            logKoaUL += 1**2
        elif solutedependencies[sltdes][1] == 4 and sltdes == 'L':
            logKoaUL += 2**2
        elif solutedependencies[sltdes][1] > 4:
            logKoaUL += 3**2
    logKoaUL = int(np.ceil((logKoaUL/4)**0.5))
    logKoaerr = (solutedependencies['V'][0] * 0.08)**2 + \
                (solutedependencies['L'][0] * 0.79)**2 * ((solutedependencies['L'][2] / solutedependencies['L'][0])**2 + (0.02 / 0.79)**2) + \
                0.03**2
    if solutedependencies['S'][0] != 0:
        logKoaerr += (solutedependencies['S'][0] * 0.69)**2 * ((solutedependencies['S'][2] / solutedependencies['S'][0])**2 + (0.05 / 0.69)**2)
    if solutedependencies['A'][0] != 0:
        logKoaerr += (solutedependencies['A'][0] * 3.56)**2 * ((solutedependencies['A'][2] / solutedependencies['A'][0])**2 + (0.04 / 3.56)**2)
    if solutedependencies['B'][0] != 0:
        logKoaerr += (solutedependencies['B'][0] * 0.73)**2 * ((solutedependencies['B'][2] / solutedependencies['B'][0])**2 + (0.04 / 0.73)**2)
    logKoaerr = logKoaerr ** 0.5
    domainnotes = [propagated_domain_notes]
    if logKoaUL <= 1:
        domainnotes.append('aggregate solute descriptor UL is in the AD')
    else:
        domainnotes.append('aggregate solute descriptor UL is out of the AD')

    return round(logKoa, round_digits), logKoaUL, round(logKoaerr, round_digits), ', '.join(domainnotes), citation, units

