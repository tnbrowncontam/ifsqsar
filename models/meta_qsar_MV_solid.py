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
chemical_inputs = {'solute min': 1, 'solute max': 1,
                   'solvent min': 0, 'solvent max': 0,
                   'component min': 0, 'component max': 0,
                   'total min': 1, 'total max': 1}
solute_dependencies_list = ['MVmlrx', 'MVmlr', 'MVmlrRings', 'MW']
solvent_dependencies_list = []
component_dependencies_list = []
propagated_domain_notes = ''
smiles_flag = 'neutrals'

stored = {}

def calculate(solutedependencies, solventdependencies, componentdependencies, solutef, solventf, componentf):
    # Sum the different fragment types defined in Kotomin and Kozolov 2006
    MV = solutedependencies[0]['MVmlrx'][0] + solutedependencies[0]['MVmlr'][0] + solutedependencies[0]['MVmlrRings'][0]
    # if there are no fragments output MV corrected to get liquid density of methane
    if MV == 0:
        MV = (solutedependencies[0]['MW'][0] / 0.4228) * (0.371/0.4228)
    return round(MV, round_digits), np.nan, round(0, round_digits), propagated_domain_notes, citation, units, endpoint

