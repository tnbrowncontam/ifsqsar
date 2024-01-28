"""Meta QSAR for logKoo (wet vs. dry octanol partitioning)"""
import numpy as np
value_names = ('logKoo',)
version = 1
endpoint = 'Log of wet octanol - dry octanol partition coefficient, used as a conversion factor'
citation = 'Solute descriptors: '\
           'Brown, T. N.; '\
           'QSPRs for Predicting Equilibrium Partitioning in Solvent-Air Systems '\
           'from the Chemical Structures of Solutes and Solvents. '\
           'J Solution Chem 2022, (https://doi.org/10.1007/s10953-022-01162-2).'\
           'PPLFER Equation: '\
           'Brown T.N., Sangion A., Arnot J.A.; '\
           'Identifying Uncertainty in Physical-Chemical Property Estimation with IFSQSAR.'\
           '2024, In Prep.'
round_digits = 2
units = 'log L[dry o]/L[wet o]'
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
    # calculate log Kow with pplfer equation from Brown 2021
    logKoo = (-1.36-(-1.57)) * solutedependencies[0]['S'][0] \
             + (-0.13-(-0.16)) * solutedependencies[0]['A'][0] \
             + (-3.49-(-4.05)) * solutedependencies[0]['B'][0] \
             + (2.41-(2.71)) * solutedependencies[0]['V'][0] \
             + (0.41-(0.41)) * solutedependencies[0]['L'][0] \
             + (0.41-(0.38))
    logKooUL = 0
    ecount = 0
    ucount = 0
    for sltdes in ['S', 'A', 'B', 'L']:
        if solutedependencies[0][sltdes][1] == 'E':
            ecount += 1
        elif solutedependencies[0][sltdes][1] == 'U':
            ucount += 1
        elif solutedependencies[0][sltdes][1] < 4:
            logKooUL += solutedependencies[0][sltdes][1]**2
        elif solutedependencies[0][sltdes][1] == 4 and sltdes != 'L':
            logKooUL += 1**2
        elif solutedependencies[0][sltdes][1] == 4 and sltdes == 'L':
            logKooUL += 2**2
        elif solutedependencies[0][sltdes][1] == 6:
            logKooUL += 3**2
    if ecount+ucount < 4:
        logKooUL = int(np.ceil((logKooUL/(4-ecount-ucount))**0.5))
    if solutedependencies[0]['S'][1] == 5 or solutedependencies[0]['A'][1] == 5 or solutedependencies[0]['B'][1] == 5 or solutedependencies[0]['L'][1] == 5:
        logKooUL = 5
    logKooerr = solutedependencies[0]['V'][0]**2 * (0.06**2 + (0.08**2 + 0.06**2)) + (0.03**2 + (0.03**2 + 0.03**2))
    if solutedependencies[0]['S'][0] != 0:
        logKooerr += (solutedependencies[0]['S'][0] * (-1.35-(-1.57)))**2 * ((solutedependencies[0]['S'][2] / solutedependencies[0]['S'][0])**2 + (0.04**2 + (0.05**2 + 0.05**2)) / (-1.35-(-1.57))**2)
    if solutedependencies[0]['A'][0] != 0:
        logKooerr += (solutedependencies[0]['A'][0] * (-0.13-(-0.16)))**2 * ((solutedependencies[0]['A'][2] / solutedependencies[0]['A'][0])**2 + (0.03**2 + (0.04**2 + 0.04**2)) / (-0.13-(-0.16))**2)
    if solutedependencies[0]['B'][0] != 0:
        logKooerr += (solutedependencies[0]['B'][0] * (-3.49-(-4.05)))**2 * ((solutedependencies[0]['B'][2] / solutedependencies[0]['B'][0])**2 + (0.03**2 + (0.04**2 + 0.04**2)) / (-3.49-(-4.05))**2)
    logKooerr = logKooerr ** 0.5
    domainnotes = [propagated_domain_notes]
    if ecount+ucount > 0:
        if ecount+ucount < 4:
            noteconcat = []
            logKooULconcat = []
            if ecount > 0:
                logKooULconcat.append('E')
                noteconcat.append('experimental')
            if ucount > 0:
                logKooULconcat.append('U')
                noteconcat.append('user')
            logKooULconcat.append(str(logKooUL))
            noteconcat.append('predicted values')
            if logKooUL <= 1:
                noteconcat.append('aggregate solute descriptor UL is in the AD')
            else:
                noteconcat.append('aggregate solute descriptor UL is out of the AD')
            logKooUL = ''.join(logKooULconcat)
            domainnotes.append(', '.join(noteconcat))
        else:
            logKooULconcat = []
            if ecount > 0:
                logKooULconcat.append('E')
            if ucount > 0:
                logKooULconcat.append('U')
            logKooUL = ''.join(logKooULconcat)
            domainnotes.append('experimental or user values, aggregate solute descriptor UL is in the AD')
    elif logKooUL <= 1:
        domainnotes.append('aggregate solute descriptor UL is in the AD')
    else:
        domainnotes.append('aggregate solute descriptor UL is out of the AD')
    if logKooUL == 'E':
        errorscale = 1
    else:
        errorscale = 1.25
    logKooerr *= errorscale

    return round(logKoo, round_digits), logKooUL, round(logKooerr, round_digits), ', '.join(domainnotes), citation, units, endpoint

