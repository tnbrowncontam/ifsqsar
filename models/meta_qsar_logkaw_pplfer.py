"""Meta QSAR for logKaw"""
import numpy as np
value_names = ('logKaw',)
version = 1
round_digits = 2
units = 'unitless'
components = {'solute': 1, 'solvent': 0}
solute_dependencies_list = ['S', 'A', 'B', 'V', 'L']
solvent_dependencies_list = []
smiles_flag = 'neutrals'


def calculate(solutedependencies, solventdependencies):
    # calculate log Koa with pplfer equation from Brown 2021
    logKaw = -2.26 * solutedependencies['S'][0] \
             - 3.72 * solutedependencies['A'][0] \
             - 4.78 * solutedependencies['B'][0] \
             + 2.19 * solutedependencies['V'][0] \
             - 0.38 * solutedependencies['L'][0] \
             + 0.64
    logKawUL = 0
    for sltdes in ['S', 'A', 'B', 'L']:
        if solutedependencies[sltdes][1] < 4:
            logKawUL += solutedependencies[sltdes][1]**2
        elif solutedependencies[sltdes][1] == 4 and sltdes != 'L':
            logKawUL += 1**2
        elif solutedependencies[sltdes][1] == 4 and sltdes == 'L':
            logKawUL += 2**2
        elif solutedependencies[sltdes][1] > 4:
            logKawUL += 3**2
    logKawUL = np.ceil((logKawUL/4)**0.5)
    logKawerr = (solutedependencies['V'][0] * 0.06)**2 + \
                (solutedependencies['L'][0] * -0.38)**2 * ((solutedependencies['L'][2] / solutedependencies['L'][0])**2 + (0.02 / -0.38)**2) + \
                0.03**2
    if solutedependencies['S'][0] != 0:
        logKawerr += (solutedependencies['S'][0] * -2.26)**2 * ((solutedependencies['S'][2] / solutedependencies['S'][0])**2 + (0.05 / -2.26)**2)
    if solutedependencies['A'][0] != 0:
        logKawerr += (solutedependencies['A'][0] * -3.72)**2 * ((solutedependencies['A'][2] / solutedependencies['A'][0])**2 + (0.04 / -3.72)**2)
    if solutedependencies['B'][0] != 0:
        logKawerr += (solutedependencies['B'][0] * -4.78)**2 * ((solutedependencies['B'][2] / solutedependencies['B'][0])**2 + (0.04 / -4.78)**2)
    logKawerr = logKawerr ** 0.5
    domainnotes = []
    if logKawUL <= 1:
        domainnotes.append('aggregate solute descriptor UL is in the AD')
    else:
        domainnotes.append('aggregate solute descriptor UL is out of the AD')

    return round(logKaw, round_digits), logKawUL, round(logKawerr, round_digits), ', '.join(domainnotes)

