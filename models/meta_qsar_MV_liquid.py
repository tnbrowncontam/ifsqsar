"""Clone of QSPR for Molar Volume of liquids (MVliquid)"""
import numpy as np
value_names = ('MVliquid',)
version = 1
endpoint = 'Molar volume of liquid'
citation = 'Kotomin, A. A.; Kozlov, A. S., '\
           'Calculation of densities of organic compounds from contributions of molecular fragments. '\
           'Russ J Appl Chem 2006, 79 (6), 957-966.'
round_digits = 2
units = 'cm^3/mol'
chemical_inputs = {'solute min': 1, 'solute max': 1,
                   'solvent min': 0, 'solvent max': 0,
                   'component min': 0, 'component max': 0,
                   'total min': 1, 'total max': 1}
solute_dependencies_list = ['MVsolid', 'MVliqcorr']
solvent_dependencies_list = []
component_dependencies_list = []
propagated_domain_notes = ''
smiles_flag = 'neutrals'

stored = {'O': (round(18.01528, round_digits), 'E', np.nan, 'experimental value used', 'well known value', units, endpoint)}

def calculate(solutedependencies, solventdependencies, componentdependencies, solutef, solventf, componentf):
    # Apply solid -> liquid correction defined in Kotomin and Kozolov 2006
    MV = solutedependencies[0]['MVsolid'][0] * solutedependencies[0]['MVliqcorr'][0]

    return round(MV, round_digits), np.nan, round(0, round_digits), propagated_domain_notes, citation, units, endpoint

