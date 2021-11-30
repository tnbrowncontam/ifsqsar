"""Meta QSAR for MVsolid"""
value_names = ('MVsolid',)
version = 1
round_digits = 2
units = 'unitless'
components = {'solute': 1, 'solvent': 0}
solute_dependencies_list = ['MVmlrx', 'MVmlr', 'MVmlrRings']
solvent_dependencies_list = []
smiles_flag = 'neutrals'


def calculate(solutedependencies, solventdependencies):
    # Sum the different fragment types defined in Kotomin and Kozolov 2006
    MV = solutedependencies['MVmlrx'][0] + solutedependencies['MVmlr'][0] + solutedependencies['MVmlrRings'][0]

    return round(MV, round_digits), 0, round(0, round_digits), ''

