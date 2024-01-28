"""Meta QSAR for logSwliquid"""
import numpy as np
value_names = ('logSwliquid',)
version = 2
endpoint = 'Log of solubility for liquids or super-cooled liquid solutes in water predicted by PPLFER at 298K'
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
units = 'log mol/L[w]'
chemical_inputs = {'solute min': 1, 'solute max': 1,
                   'solvent min': 0, 'solvent max': 0,
                   'component min': 0, 'component max': 0,
                   'total min': 1, 'total max': 1}
solute_dependencies_list = ['S', 'A', 'B', 'V', 'L', 'MVliquid', 'state']
solvent_dependencies_list = []
component_dependencies_list = []
propagated_domain_notes = ''
smiles_flag = 'neutrals'

stored = {}

def calculate(solutedependencies, solventdependencies, componentdependencies, solutef, solventf, componentf):
    domainnotes = [propagated_domain_notes]
    # calculate AB and ABerr
    AB = (solutedependencies[0]['A'][0] * solutedependencies[0]['B'][0])**0.5
    if solutedependencies[0]['A'][0] > 0 and solutedependencies[0]['B'][0] > 0:
        ABerr = solutedependencies[0]['A'][0] * solutedependencies[0]['B'][0] * ((solutedependencies[0]['A'][2]/solutedependencies[0]['A'][0])**2 + (solutedependencies[0]['B'][2]/solutedependencies[0]['B'][0])**2)**0.5
        ABerr = AB * 0.5 * ABerr / (solutedependencies[0]['A'][0] * solutedependencies[0]['B'][0])
    else:
        ABerr = 0
    # calculate logSw with pplfer equation from Brown et al. 2023
    logSw = 0.71 * solutedependencies[0]['S'][0] \
            + 2.80 * solutedependencies[0]['A'][0] \
            + 4.15 * solutedependencies[0]['B'][0] \
            - 1.60 * AB \
            - 3.49 * solutedependencies[0]['V'][0] \
            - 0.13 * solutedependencies[0]['L'][0] \
            + 0.18
    # calculate UL and err
    logSwUL = 0
    ecount = 0
    ucount = 0
    for sltdes in ['S', 'A', 'B', 'L']:
        if solutedependencies[0][sltdes][1] == 'E':
            ecount += 1
        elif solutedependencies[0][sltdes][1] == 'U':
            ucount += 1
        elif solutedependencies[0][sltdes][1] < 4:
            logSwUL += solutedependencies[0][sltdes][1]**2
        elif solutedependencies[0][sltdes][1] == 4 and sltdes != 'L':
            logSwUL += 1**2
        elif solutedependencies[0][sltdes][1] == 4 and sltdes == 'L':
            logSwUL += 2**2
        elif solutedependencies[0][sltdes][1] == 6:
            logSwUL += 3**2
    if ecount+ucount < 4:
        logSwUL = int(np.ceil((logSwUL/(4-ecount-ucount))**0.5))
    if solutedependencies[0]['S'][1] == 5 or solutedependencies[0]['A'][1] == 5 or solutedependencies[0]['B'][1] == 5 or solutedependencies[0]['L'][1] == 5:
        logSwUL = 5
    logSwerr = (solutedependencies[0]['V'][0] * 0.17)**2 + \
                (solutedependencies[0]['L'][0] * -0.13)**2 * ((solutedependencies[0]['L'][2] / solutedependencies[0]['L'][0])**2 + (0.04 / -0.13)**2) + \
                0.07**2
    if solutedependencies[0]['S'][0] != 0:
        logSwerr += (solutedependencies[0]['S'][0] * 0.71)**2 * ((solutedependencies[0]['S'][2] / solutedependencies[0]['S'][0])**2 + (0.11 / 0.71)**2)
    if solutedependencies[0]['A'][0] != 0:
        logSwerr += (solutedependencies[0]['A'][0] * 2.80)**2 * ((solutedependencies[0]['A'][2] / solutedependencies[0]['A'][0])**2 + (0.23 / 2.80)**2)
    if solutedependencies[0]['B'][0] != 0:
        logSwerr += (solutedependencies[0]['B'][0] * 4.15)**2 * ((solutedependencies[0]['B'][2] / solutedependencies[0]['B'][0])**2 + (0.13 / 4.15)**2)
    if AB != 0:
        logSwerr += (AB * -1.60)**2 * ((ABerr / AB)**2 + (0.27 / -1.60)**2)
    logSwerr = logSwerr ** 0.5
    if ecount+ucount > 0:
        if ecount+ucount < 4:
            noteconcat = []
            logSwULconcat = []
            if ecount > 0:
                logSwULconcat.append('E')
                noteconcat.append('experimental')
            if ucount > 0:
                logSwULconcat.append('U')
                noteconcat.append('user')
            logSwULconcat.append(str(logSwUL))
            noteconcat.append('predicted values')
            if logSwUL <= 1:
                noteconcat.append('aggregate solute descriptor UL is in the AD')
            else:
                noteconcat.append('aggregate solute descriptor UL is out of the AD')
            logSwUL = ''.join(logSwULconcat)
            domainnotes.append(', '.join(noteconcat))
        else:
            logSwULconcat = []
            if ecount > 0:
                logSwULconcat.append('E')
            if ucount > 0:
                logSwULconcat.append('U')
            logSwUL = ''.join(logSwULconcat)
            domainnotes.append('experimental or user values, aggregate solute descriptor UL is in the AD')
    elif logSwUL <= 1:
        domainnotes.append('aggregate solute descriptor UL is in the AD')
    else:
        domainnotes.append('aggregate solute descriptor UL is out of the AD')
    if logSwUL == 'E':
        errorscale = 1
    else:
        errorscale = 1.25
        if solutedependencies[0]['state'][3] in ['likely solid', 'maybe solid']:
            errorscale *= 1.25
    logSwerr *= errorscale
    # cap Sw at the inverse of solute MV
    logSwmax = np.log10(1000 / solutedependencies[0]['MVliquid'][0])
    if logSw > logSwmax:
        domainnotes = ['Predicted water solubility ({}) capped at inverse of molar volume, UL set to 6; original aggregate UL: {}'.format(round(logSw, round_digits), logSwUL)] + domainnotes
        logSw = logSwmax
        logSwUL = 6
    return round(logSw, round_digits), logSwUL, round(logSwerr, round_digits), '; '.join(domainnotes), citation, units, endpoint

