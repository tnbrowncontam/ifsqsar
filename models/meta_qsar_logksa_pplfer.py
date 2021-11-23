"""Meta QSAR for logKsa"""
value_names = ('logKsa',)
version = 1
round_digits = 2
units = 'unitless'
components = {'solute': 1, 'solvent': 1}
solute_dependencies_list = ['E', 'S', 'A', 'B', 'V', 'L']
solvent_dependencies_list = ['s', 'a', 'b', 'v', 'l', 'c']
smiles_flag = 'neutrals'


def calculate(solutedependencies, solventdependencies):
    # calculate log Kow with pplfer equation from Brown 2021
    logKsa = 0

    return round(logKsa, round_digits), 0, round(0, round_digits), ''

