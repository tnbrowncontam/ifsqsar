"""Meta QSAR for logKaw"""
import numpy as np
value_names = ('logKaw',)
version = 1
endpoint = 'Log of air-water partition coefficient (Henry\'s Law Constant)'
citation = 'Solute descriptors: '\
           'Brown, T. N.; '\
           'QSPRs for Predicting Equilibrium Partitioning in Solvent-Air Systems '\
           'from the Chemical Structures of Solutes and Solvents. '\
           'J Solution Chem 2022, (https://doi.org/10.1007/s10953-022-01162-2).'\
           'PPLFER Equation: '\
           'Brown, T. N., '\
           'Empirical Regressions between System Parameters and Solute Descriptors of Polyparameter Linear Free '\
           'Energy Relationships (PPLFERs) for Predicting Solvent-Air Partitioning. '\
           'Fluid Phase Equilibria 2021, 113035.'
round_digits = 2
units = 'log L[w]/L[a]'
chemical_inputs = {'solute min': 1, 'solute max': 1,
                   'solvent min': 0, 'solvent max': 0,
                   'component min': 0, 'component max': 0,
                   'total min': 1, 'total max': 1}
solute_dependencies_list = ['S', 'A', 'B', 'V', 'L']
solvent_dependencies_list = []
component_dependencies_list = []
propagated_domain_notes = ''
smiles_flag = 'neutrals'

stored = {}

def calculate(solutedependencies, solventdependencies, componentdependencies, solutef, solventf, componentf):
    # calculate log Koa with pplfer equation from Brown 2021
    logKaw = -2.26 * solutedependencies[0]['S'][0] \
             - 3.72 * solutedependencies[0]['A'][0] \
             - 4.78 * solutedependencies[0]['B'][0] \
             + 2.19 * solutedependencies[0]['V'][0] \
             - 0.38 * solutedependencies[0]['L'][0] \
             + 0.64
    logKawUL = 0
    ecount = 0
    ucount = 0
    for sltdes in ['S', 'A', 'B', 'L']:
        if solutedependencies[0][sltdes][1] == 'E':
            ecount += 1
        elif solutedependencies[0][sltdes][1] == 'U':
            ucount += 1
        elif solutedependencies[0][sltdes][1] < 4:
            logKawUL += solutedependencies[0][sltdes][1]**2
        elif solutedependencies[0][sltdes][1] == 4 and sltdes != 'L':
            logKawUL += 1**2
        elif solutedependencies[0][sltdes][1] == 4 and sltdes == 'L':
            logKawUL += 2**2
        elif solutedependencies[0][sltdes][1] == 6:
            logKawUL += 3**2
    if ecount+ucount < 4:
        logKawUL = int(np.ceil((logKawUL/(4-ecount-ucount))**0.5))
    if solutedependencies[0]['S'][1] == 5 or solutedependencies[0]['A'][1] == 5 or solutedependencies[0]['B'][1] == 5 or solutedependencies[0]['L'][1] == 5:
        logKawUL = 5
    logKawerr = (solutedependencies[0]['V'][0] * 0.06)**2 + \
                (solutedependencies[0]['L'][0] * -0.38)**2 * ((solutedependencies[0]['L'][2] / solutedependencies[0]['L'][0])**2 + (0.02 / -0.38)**2) + \
                0.03**2
    if solutedependencies[0]['S'][0] != 0:
        logKawerr += (solutedependencies[0]['S'][0] * -2.26)**2 * ((solutedependencies[0]['S'][2] / solutedependencies[0]['S'][0])**2 + (0.05 / -2.26)**2)
    if solutedependencies[0]['A'][0] != 0:
        logKawerr += (solutedependencies[0]['A'][0] * -3.72)**2 * ((solutedependencies[0]['A'][2] / solutedependencies[0]['A'][0])**2 + (0.04 / -3.72)**2)
    if solutedependencies[0]['B'][0] != 0:
        logKawerr += (solutedependencies[0]['B'][0] * -4.78)**2 * ((solutedependencies[0]['B'][2] / solutedependencies[0]['B'][0])**2 + (0.04 / -4.78)**2)
    logKawerr = logKawerr ** 0.5
    domainnotes = [propagated_domain_notes]
    if ecount+ucount > 0:
        if ecount+ucount < 4:
            noteconcat = []
            logKawULconcat = []
            if ecount > 0:
                logKawULconcat.append('E')
                noteconcat.append('experimental')
            if ucount > 0:
                logKawULconcat.append('U')
                noteconcat.append('user')
            logKawULconcat.append(str(logKawUL))
            noteconcat.append('predicted values')
            if logKawUL <= 1:
                noteconcat.append('aggregate solute descriptor UL is in the AD')
            else:
                noteconcat.append('aggregate solute descriptor UL is out of the AD')
            logKawUL = ''.join(logKawULconcat)
            domainnotes.append(', '.join(noteconcat))
        else:
            logKawULconcat = []
            if ecount > 0:
                logKawULconcat.append('E')
            if ucount > 0:
                logKawULconcat.append('U')
            logKawUL = ''.join(logKawULconcat)
            domainnotes.append('experimental or user values, aggregate solute descriptor UL is in the AD')
    elif logKawUL <= 1:
        domainnotes.append('aggregate solute descriptor UL is in the AD')
    else:
        domainnotes.append('aggregate solute descriptor UL is out of the AD')
    if logKawUL == 'E':
        errorscale = 1
    else:
        errorscale = 1.25
    logKawerr *= errorscale

    return round(logKaw, round_digits), logKawUL, round(logKawerr, round_digits), '; '.join(domainnotes), citation, units, endpoint

