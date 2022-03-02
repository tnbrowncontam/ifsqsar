"""Meta QSAR for logKow (dry solvent)"""
import numpy as np
value_names = ('logKow',)
version = 1
endpoint = 'Log of octanol-water partition coefficient (log P)'
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
units = 'log L[w]/L[o]'
components = {'solute': 1, 'solvent': 0}
solute_dependencies_list = ['S', 'A', 'B', 'V', 'L']
solvent_dependencies_list = []
propagated_domain_notes = ''
smiles_flag = 'neutrals'


def calculate(solutedependencies, solventdependencies):
    # calculate log Kow with pplfer equation from Brown 2021
    logKow = (0.69-2.26) * solutedependencies['S'][0] \
             + (3.56-3.72) * solutedependencies['A'][0] \
             + (0.73-4.78) * solutedependencies['B'][0] \
             + (0.52-(-2.19)) * solutedependencies['V'][0] \
             + (0.79-0.38) * solutedependencies['L'][0] \
             + (-0.26-(-0.64))
    logKowUL = 0
    for sltdes in ['S', 'A', 'B', 'L']:
        if solutedependencies[sltdes][1] < 4:
            logKowUL += solutedependencies[sltdes][1]**2
        elif solutedependencies[sltdes][1] == 4 and sltdes != 'L':
            logKowUL += 1**2
        elif solutedependencies[sltdes][1] == 4 and sltdes == 'L':
            logKowUL += 2**2
        elif solutedependencies[sltdes][1] > 4:
            logKowUL += 3**2
    logKowUL = int(np.ceil((logKowUL/4)**0.5))
    logKowerr = solutedependencies['V'][0]**2 * (0.08**2 + 0.06**2) + \
                (solutedependencies['L'][0] * (0.79-0.38))**2 * ((solutedependencies['L'][2] / solutedependencies['L'][0])**2 + (0.02**2 + 0.02**2) / (0.79-0.38)**2) + \
                (0.03**2 + 0.03**2)
    if solutedependencies['S'][0] != 0:
        logKowerr += (solutedependencies['S'][0] * (0.69-2.26))**2 * ((solutedependencies['S'][2] / solutedependencies['S'][0])**2 + (0.05**2 + 0.05**2) / (0.69-2.26)**2)
    if solutedependencies['A'][0] != 0:
        logKowerr += (solutedependencies['A'][0] * (3.56-3.72))**2 * ((solutedependencies['A'][2] / solutedependencies['A'][0])**2 + (0.04**2 + 0.04**2) / (3.56-3.72)**2)
    if solutedependencies['B'][0] != 0:
        logKowerr += (solutedependencies['B'][0] * (0.73-4.78))**2 * ((solutedependencies['B'][2] / solutedependencies['B'][0])**2 + (0.04**2 + 0.04**2) / (0.73-4.78)**2)
    logKowerr = logKowerr ** 0.5
    domainnotes = [propagated_domain_notes]
    if logKowUL <= 1:
        domainnotes.append('aggregate solute descriptor UL is in the AD')
    else:
        domainnotes.append('aggregate solute descriptor UL is out of the AD')

    return round(logKow, round_digits), logKowUL, round(logKowerr, round_digits), ', '.join(domainnotes), citation, units, endpoint

