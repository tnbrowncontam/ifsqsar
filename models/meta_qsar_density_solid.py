"""Clone of QSPR for density of solids (densitysolid)"""
value_names = ('densitysolid',)
version = 1
citation = 'Kotomin, A. A.; Kozlov, A. S., '\
           'Calculation of densities of organic compounds from contributions of molecular fragments. '\
           'Russ J Appl Chem 2006, 79 (6), 957-966.'
round_digits = 3
units = 'g/cm^3'
components = {'solute': 1, 'solvent': 0}
solute_dependencies_list = ['MVsolid', 'MW']
solvent_dependencies_list = []
smiles_flag = 'neutrals'


def calculate(solutedependencies, solventdependencies):
    density = solutedependencies['MW'][0] / solutedependencies['MVsolid'][0]

    return round(density, round_digits), 0, round(0, round_digits), '', citation

