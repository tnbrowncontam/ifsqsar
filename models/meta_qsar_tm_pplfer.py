"""Meta QSAR for tmpplfer"""
import numpy as np
value_names = ('tmpplfer',)
version = 1
round_digits = 2
units = 'K'
components = {'solute': 1, 'solvent': 0}
solute_dependencies_list = ['E', 'S', 'A', 'B', 'L']
solvent_dependencies_list = []
smiles_flag = 'neutrals'

def calculate(solutedependencies, solventdependencies):
    domainnotes = []
    # calculate MP
    MP = 53.6 * solutedependencies['E'][0] + \
         27.8 * solutedependencies['S'][0] + \
         107.9 * solutedependencies['A'][0] + \
         -19.2 * solutedependencies['B'][0] + \
         7.1 * solutedependencies['L'][0] + \
         -90.6 + \
         273.15
    # calculate aggregate UL
    MPUL = 0
    for sltdes in ['E', 'S', 'A', 'B', 'L']:
        if solutedependencies[sltdes][1] < 4:
            MPUL += solutedependencies[sltdes][1]**2
        elif solutedependencies[sltdes][1] == 4 and sltdes != 'L':
            MPUL += 1**2
        elif solutedependencies[sltdes][1] == 4 and sltdes == 'L':
            MPUL += 2**2
        elif solutedependencies[sltdes][1] > 4:
            MPUL += 3**2
    MPUL = np.ceil((MPUL/5)**0.5)
    # calculate error
    MPerr = (solutedependencies['L'][0] * 7.1)**2 * ((solutedependencies['L'][2] / solutedependencies['L'][0])**2 + (0.7 / 7.1)**2) + \
            3.0**2
    if solutedependencies['E'][0] != 0:
        MPerr += (solutedependencies['E'][0] * 53.6)**2 * ((solutedependencies['E'][2] / solutedependencies['E'][0])**2 + (3.6 / 53.6)**2)
    if solutedependencies['S'][0] != 0:
        MPerr += (solutedependencies['S'][0] * 27.8)**2 * ((solutedependencies['S'][2] / solutedependencies['S'][0])**2 + (4.1 / 27.8)**2)
    if solutedependencies['A'][0] != 0:
        MPerr += (solutedependencies['A'][0] * 107.9)**2 * ((solutedependencies['A'][2] / solutedependencies['A'][0])**2 + (4.9 / 107.9)**2)
    if solutedependencies['B'][0] != 0:
        MPerr += (solutedependencies['B'][0] * -19.2)**2 * ((solutedependencies['B'][2] / solutedependencies['B'][0])**2 + (4.4 / -19.2)**2)
    MPerr = MPerr**0.5

    return round(MP, round_digits), int(MPUL), round(MPerr, round_digits), ', '.join(domainnotes)