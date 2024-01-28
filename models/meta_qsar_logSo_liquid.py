"""Meta QSAR for logSoliquid (dry solvent)"""
import numpy as np
value_names = ('logSoliquid',)
version = 1
endpoint = 'Log of solubility in dry octanol for liquid or super-cooled liquid solute'
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
units = 'log mol/L'
chemical_inputs = {'solute min': 1, 'solute max': 1,
                   'solvent min': 0, 'solvent max': 0,
                   'component min': 0, 'component max': 0,
                   'total min': 1, 'total max': 1}
solute_dependencies_list = ['S', 'A', 'B', 'V', 'L', 'MVliquid']
solvent_dependencies_list = []
component_dependencies_list = []
propagated_domain_notes = ''
smiles_flag = 'neutrals'

stored = {}

def calculate(solutedependencies, solventdependencies, componentdependencies, solutef, solventf, componentf):
    domainnotes = [propagated_domain_notes]
    phaseerrorscaling = 1
    # calculate AB and ABerr
    AB = (solutedependencies[0]['A'][0] * solutedependencies[0]['B'][0])**0.5
    if solutedependencies[0]['A'][0] > 0 and solutedependencies[0]['B'][0] > 0:
        ABerr = solutedependencies[0]['A'][0] * solutedependencies[0]['B'][0] * ((solutedependencies[0]['A'][2]/solutedependencies[0]['A'][0])**2 + (solutedependencies[0]['B'][2]/solutedependencies[0]['B'][0])**2)**0.5
        ABerr = AB * 0.5 * ABerr / (solutedependencies[0]['A'][0] * solutedependencies[0]['B'][0])
    else:
        ABerr = 0
    # calculate logSo with pplfer equation from Brown et al. 2023
    logSo = - 0.86 * solutedependencies[0]['S'][0] \
            + 2.64 * solutedependencies[0]['A'][0] \
            + 0.10 * solutedependencies[0]['B'][0] \
            - 1.60 * AB \
            - 0.78 * solutedependencies[0]['V'][0] \
            + 0.28 * solutedependencies[0]['L'][0] \
            + 0.56
    # calculate UL and err
    logSoUL = 0
    ecount = 0
    ucount = 0
    for sltdes in ['S', 'A', 'B', 'L']:
        if solutedependencies[0][sltdes][1] == 'E':
            ecount += 1
        elif solutedependencies[0][sltdes][1] == 'U':
            ucount += 1
        elif solutedependencies[0][sltdes][1] < 4:
            logSoUL += solutedependencies[0][sltdes][1]**2
        elif solutedependencies[0][sltdes][1] == 4 and sltdes != 'L':
            logSoUL += 1**2
        elif solutedependencies[0][sltdes][1] == 4 and sltdes == 'L':
            logSoUL += 2**2
        elif solutedependencies[0][sltdes][1] == 6:
            logSoUL += 3**2
    if ecount+ucount < 4:
        logSoUL = int(np.ceil((logSoUL/(4-ecount-ucount))**0.5))
    if solutedependencies[0]['S'][1] == 5 or solutedependencies[0]['A'][1] == 5 or solutedependencies[0]['B'][1] == 5 or solutedependencies[0]['L'][1] == 5:
        logSoUL = 5
    logSoerr = (solutedependencies[0]['V'][0] * 0.20)**2 + \
                (solutedependencies[0]['L'][0] * 0.28)**2 * ((solutedependencies[0]['L'][2] / solutedependencies[0]['L'][0])**2 + (0.05 / 0.28)**2) + \
                0.08**2
    if solutedependencies[0]['S'][0] != 0:
        logSoerr += (solutedependencies[0]['S'][0] * -0.86)**2 * ((solutedependencies[0]['S'][2] / solutedependencies[0]['S'][0])**2 + (0.13 / -0.86)**2)
    if solutedependencies[0]['A'][0] != 0:
        logSoerr += (solutedependencies[0]['A'][0] * 2.64)**2 * ((solutedependencies[0]['A'][2] / solutedependencies[0]['A'][0])**2 + (0.23 / 2.64)**2)
    if solutedependencies[0]['B'][0] != 0:
        logSoerr += (solutedependencies[0]['B'][0] * 0.10)**2 * ((solutedependencies[0]['B'][2] / solutedependencies[0]['B'][0])**2 + (0.14 / 0.10)**2)
    if AB != 0:
        logSoerr += (AB * -1.60)**2 * ((ABerr / AB)**2 + (0.27 / -1.60)**2)
    logSoerr = logSoerr ** 0.5
    if ecount+ucount > 0:
        if ecount+ucount < 4:
            noteconcat = []
            logSoULconcat = []
            if ecount > 0:
                logSoULconcat.append('E')
                noteconcat.append('experimental')
            if ucount > 0:
                logSoULconcat.append('U')
                noteconcat.append('user')
            logSoULconcat.append(str(logSoUL))
            noteconcat.append('predicted values')
            if logSoUL <= 1:
                noteconcat.append('aggregate solute descriptor UL is in the AD')
            else:
                noteconcat.append('aggregate solute descriptor UL is out of the AD')
            logSoUL = ''.join(logSoULconcat)
            domainnotes.append(', '.join(noteconcat))
        else:
            logSoULconcat = []
            if ecount > 0:
                logSoULconcat.append('E')
            if ucount > 0:
                logSoULconcat.append('U')
            logSoUL = ''.join(logSoULconcat)
            domainnotes.append('experimental or user values, aggregate solute descriptor UL is in the AD')
    elif logSoUL <= 1:
        domainnotes.append('aggregate solute descriptor UL is in the AD')
    else:
        domainnotes.append('aggregate solute descriptor UL is out of the AD')
    if logSoUL == 'E':
        errorscale = 1
    else:
        errorscale = 1.25
    logSoerr *= errorscale
    # cap So at the inverse of solute MV
    logSomax = np.log10(1000 / solutedependencies[0]['MVliquid'][0])
    if logSo > logSomax:
        domainnotes = ['Predicted dry octanol solubility ({}) capped at inverse of molar volume, UL set to 6; original aggregate UL: {}'.format(round(logSo, round_digits), logSoUL)] + domainnotes
        logSo = logSomax
        logSoUL = 6
    return round(logSo, round_digits), logSoUL, round(phaseerrorscaling * logSoerr, round_digits), '; '.join(domainnotes), citation, units, endpoint

