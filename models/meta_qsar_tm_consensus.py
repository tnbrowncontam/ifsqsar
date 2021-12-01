"""Meta QSAR for tmconsensus"""
import numpy as np
value_names = ('tmconsensus',)
version = 1
round_digits = 2
units = 'K'
components = {'solute': 1, 'solvent': 0}
solute_dependencies_list = ['tm', 'tmpplfer']
solvent_dependencies_list = []
smiles_flag = 'neutrals'

def calculate(solutedependencies, solventdependencies):
    domainnotes = []
    MP = (solutedependencies['tm'][0] + solutedependencies['tmpplfer'][0]) / 2
    MPUL = np.ceil(((solutedependencies['tm'][1]**2 + solutedependencies['tmpplfer'][1]**2) / 2)**0.5)
    MPerr = ((0.5 * solutedependencies['tm'][2])**2 + (0.5 * solutedependencies['tmpplfer'][2])**2)**0.5
    domainnotes.append('tmqsar notes:')
    domainnotes.append(solutedependencies['tm'][3])
    domainnotes.append('tmpplfer notes:')
    domainnotes.append(solutedependencies['tmpplfer'][3])
    return round(MP, round_digits), int(MPUL), round(MPerr, round_digits), ', '.join(domainnotes)

