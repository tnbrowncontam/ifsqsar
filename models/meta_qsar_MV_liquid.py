"""Clone of QSPR for Molar Volume of liquids (MVliquid)"""
value_names = ('MVliquid',)
version = 1
citation = 'Kotomin, A. A.; Kozlov, A. S., '\
           'Calculation of densities of organic compounds from contributions of molecular fragments. '\
           'Russ J Appl Chem 2006, 79 (6), 957-966.'
round_digits = 2
units = 'cm^3/mol'
components = {'solute': 1, 'solvent': 0}
solute_dependencies_list = ['MVsolid', 'MVliqcorr']
solvent_dependencies_list = []
smiles_flag = 'neutrals'


def calculate(solutedependencies, solventdependencies):
    # Apply solid -> liquid correction defined in Kotomin and Kozolov 2006
    MV = solutedependencies['MVsolid'][0] * solutedependencies['MVliqcorr'][0]

    return round(MV, round_digits), 0, round(0, round_digits), '', citation

