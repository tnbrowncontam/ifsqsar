"""Meta QSAR for tbpplfer"""
import numpy as np
value_names = ('tbpplfer',)
version = 1
round_digits = 2
units = 'K'
components = {'solute': 1, 'solvent': 0}
solute_dependencies_list = ['E', 'S', 'A', 'B', 'V', 'L']
solvent_dependencies_list = []
smiles_flag = 'neutrals'

def calculate(solutedependencies, solventdependencies):
    domainnotes = []
    # calculate BP
    BP = 13.0 * solutedependencies['E'][0] + \
         43.8 * solutedependencies['S'][0] + \
         59.8 * solutedependencies['A'][0] + \
         18.2 * solutedependencies['B'][0] + \
         26.9 * solutedependencies['V'][0] + \
         29.2 * solutedependencies['L'][0] + \
         -33.0 + \
         273.15
    # calculate aggregate UL
    BPUL = 0
    for sltdes in ['E', 'S', 'A', 'B', 'L']:
        if solutedependencies[sltdes][1] < 4:
            BPUL += solutedependencies[sltdes][1]**2
        elif solutedependencies[sltdes][1] == 4 and sltdes != 'L':
            BPUL += 1**2
        elif solutedependencies[sltdes][1] == 4 and sltdes == 'L':
            BPUL += 2**2
        elif solutedependencies[sltdes][1] > 4:
            BPUL += 3**2
    BPUL = np.ceil((BPUL/5)**0.5)
    # calculate error
    BPerr = (solutedependencies['V'][0] * 8.2)**2 + \
            (solutedependencies['L'][0] * 29.2)**2 * ((solutedependencies['L'][2] / solutedependencies['L'][0])**2 + (2.4 / 29.2)**2) + \
            2.8**2
    if solutedependencies['E'][0] != 0:
        BPerr += (solutedependencies['E'][0] * 13.0)**2 * ((solutedependencies['E'][2] / solutedependencies['E'][0])**2 + (3.6 / 13.0)**2)
    if solutedependencies['S'][0] != 0:
        BPerr += (solutedependencies['S'][0] * 43.8)**2 * ((solutedependencies['S'][2] / solutedependencies['S'][0])**2 + (3.3 / 43.8)**2)
    if solutedependencies['A'][0] != 0:
        BPerr += (solutedependencies['A'][0] * 59.8)**2 * ((solutedependencies['A'][2] / solutedependencies['A'][0])**2 + (3.3 / 59.8)**2)
    if solutedependencies['B'][0] != 0:
        BPerr += (solutedependencies['B'][0] * 18.2)**2 * ((solutedependencies['B'][2] / solutedependencies['B'][0])**2 + (2.8 / 18.2)**2)
    BPerr = BPerr**0.5

    return round(BP, round_digits), int(BPUL), round(BPerr, round_digits), ', '.join(domainnotes)