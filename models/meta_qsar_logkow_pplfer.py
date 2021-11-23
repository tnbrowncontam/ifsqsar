"""Meta QSAR for logKow"""
value_names = ('logKow',)
version = 1
round_digits = 2
units = 'unitless'
components = {'solute': 1, 'solvent': 0}
solute_dependencies_list = ['S', 'A', 'B', 'V', 'L']
solvent_dependencies_list = []
smiles_flag = 'neutrals'


def calculate(solutedependencies, solventdependencies):
    # calculate log Kow with pplfer equation from Brown 2021
    logKow = -1.57 * solutedependencies['S'][0] \
             - 0.16 * solutedependencies['A'][0] \
             - 4.05 * solutedependencies['B'][0] \
             + 2.71 * solutedependencies['V'][0] \
             + 0.41 * solutedependencies['L'][0] \
             + 0.38

    return round(logKow, round_digits), 0, round(0, round_digits), ''

