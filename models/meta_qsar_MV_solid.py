"""Clone of QSPR for Molar Volume of solids (MVsolid)"""
import numpy as np
value_names = ('MVsolid',)
version = 1
endpoint = 'Molar volume of solid'
citation = 'Kotomin, A. A.; Kozlov, A. S., '\
           'Calculation of densities of organic compounds from contributions of molecular fragments. '\
           'Russ J Appl Chem 2006, 79 (6), 957-966.'
round_digits = 2
units = 'cm^3/mol'
components = {'solute': 1, 'solvent': 0}
solute_dependencies_list = ['MVmlrx', 'MVmlr', 'MVmlrRings']
solvent_dependencies_list = []
propagated_domain_notes = ''
smiles_flag = 'neutrals'


def calculate(solutedependencies, solventdependencies):
    # Sum the different fragment types defined in Kotomin and Kozolov 2006
    MV = solutedependencies['MVmlrx'][0] + solutedependencies['MVmlr'][0] + solutedependencies['MVmlrRings'][0]

    return round(MV, round_digits), np.nan, round(0, round_digits), propagated_domain_notes, citation, units, endpoint

