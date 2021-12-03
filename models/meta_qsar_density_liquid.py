"""Meta QSAR for density liquid"""
value_names = ('densityliquid',)
version = 1
round_digits = 3
units = 'g/cm^3'
components = {'solute': 1, 'solvent': 0}
solute_dependencies_list = ['MVliquid', 'MW']
solvent_dependencies_list = []
smiles_flag = 'neutrals'


def calculate(solutedependencies, solventdependencies):
    density = solutedependencies['MW'][0] / solutedependencies['MVliquid'][0]

    return round(density, round_digits), 0, round(0, round_digits), ''

