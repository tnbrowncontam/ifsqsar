"""Meta QSAR for logKoa (dry solvent)"""
import numpy as np
value_names = ('logKoa',)
version = 1
endpoint = 'Log of octanol-air partition coefficient'
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
units = 'log L[a]/L[o]'
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
    logKoa = 0.69 * solutedependencies[0]['S'][0] \
             + 3.56 * solutedependencies[0]['A'][0] \
             + 0.73 * solutedependencies[0]['B'][0] \
             + 0.52 * solutedependencies[0]['V'][0] \
             + 0.79 * solutedependencies[0]['L'][0] \
             -0.26
    logKoaUL = 0
    ecount = 0
    ucount = 0
    for sltdes in ['S', 'A', 'B', 'L']:
        if solutedependencies[0][sltdes][1] == 'E':
            ecount += 1
        elif solutedependencies[0][sltdes][1] == 'U':
            ucount += 1
        elif solutedependencies[0][sltdes][1] < 4:
            logKoaUL += solutedependencies[0][sltdes][1]**2
        elif solutedependencies[0][sltdes][1] == 4 and sltdes != 'L':
            logKoaUL += 1**2
        elif solutedependencies[0][sltdes][1] == 4 and sltdes == 'L':
            logKoaUL += 2**2
        elif solutedependencies[0][sltdes][1] == 6:
            logKoaUL += 3**2
    if ecount+ucount < 4:
        logKoaUL = int(np.ceil((logKoaUL/(4-ecount-ucount))**0.5))
    if solutedependencies[0]['S'][1] == 5 or solutedependencies[0]['A'][1] == 5 or solutedependencies[0]['B'][1] == 5 or solutedependencies[0]['L'][1] == 5:
        logKoaUL = 5
    logKoaerr = (solutedependencies[0]['V'][0] * 0.08)**2 + \
                (solutedependencies[0]['L'][0] * 0.79)**2 * ((solutedependencies[0]['L'][2] / solutedependencies[0]['L'][0])**2 + (0.02 / 0.79)**2) + \
                0.03**2
    if solutedependencies[0]['S'][0] != 0:
        logKoaerr += (solutedependencies[0]['S'][0] * 0.69)**2 * ((solutedependencies[0]['S'][2] / solutedependencies[0]['S'][0])**2 + (0.05 / 0.69)**2)
    if solutedependencies[0]['A'][0] != 0:
        logKoaerr += (solutedependencies[0]['A'][0] * 3.56)**2 * ((solutedependencies[0]['A'][2] / solutedependencies[0]['A'][0])**2 + (0.04 / 3.56)**2)
    if solutedependencies[0]['B'][0] != 0:
        logKoaerr += (solutedependencies[0]['B'][0] * 0.73)**2 * ((solutedependencies[0]['B'][2] / solutedependencies[0]['B'][0])**2 + (0.04 / 0.73)**2)
    logKoaerr = logKoaerr ** 0.5
    domainnotes = [propagated_domain_notes]
    if ecount+ucount > 0:
        if ecount+ucount < 4:
            noteconcat = []
            logKoaULconcat = []
            if ecount > 0:
                logKoaULconcat.append('E')
                noteconcat.append('experimental')
            if ucount > 0:
                logKoaULconcat.append('U')
                noteconcat.append('user')
            logKoaULconcat.append(str(logKoaUL))
            noteconcat.append('predicted values')
            if logKoaUL <= 1:
                noteconcat.append('aggregate solute descriptor UL is in the AD')
            else:
                noteconcat.append('aggregate solute descriptor UL is out of the AD')
            logKoaUL = ''.join(logKoaULconcat)
            domainnotes.append(', '.join(noteconcat))
        else:
            logKoaULconcat = []
            if ecount > 0:
                logKoaULconcat.append('E')
            if ucount > 0:
                logKoaULconcat.append('U')
            logKoaUL = ''.join(logKoaULconcat)
            domainnotes.append('experimental or user values, aggregate solute descriptor UL is in the AD')
    elif logKoaUL <= 1:
        domainnotes.append('aggregate solute descriptor UL is in the AD')
    else:
        domainnotes.append('aggregate solute descriptor UL is out of the AD')
    if logKoaUL == 'E':
        errorscale = 1
    else:
        errorscale = 1.25
    logKoaerr *= errorscale

    return round(logKoa, round_digits), logKoaUL, round(logKoaerr, round_digits), ', '.join(domainnotes), citation, units, endpoint

