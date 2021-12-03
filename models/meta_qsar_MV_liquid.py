"""Meta QSAR for MVliquid"""
value_names = ('MVliquid',)
version = 1
round_digits = 2
units = 'cm^3/mol'
components = {'solute': 1, 'solvent': 0}
solute_dependencies_list = ['MVsolid', 'MVliqcorr']
solvent_dependencies_list = []
smiles_flag = 'neutrals'


def calculate(solutedependencies, solventdependencies):
    # Apply solid -> liquid correction defined in Kotomin and Kozolov 2006
    MV = solutedependencies['MVsolid'][0] * solutedependencies['MVliqcorr'][0]

    return round(MV, round_digits), 0, round(0, round_digits), ''

