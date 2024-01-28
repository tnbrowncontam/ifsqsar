"""Clone of QSPR for density of solids (densitysolid)"""
import numpy as np
value_names = ('densitysolid',)
version = 1
endpoint = 'Density of pure solid'
citation = 'Brown T.N., Sangion A., Arnot J.A.; '\
           'Identifying Uncertainty in Physical-Chemical Property Estimation with IFSQSAR.'\
           '2024, In Prep.'
round_digits = 3
units = 'g/cm^3'
chemical_inputs = {'solute min': 1, 'solute max': 1,
                   'solvent min': 0, 'solvent max': 0,
                   'component min': 0, 'component max': 0,
                   'total min': 1, 'total max': 1}
solute_dependencies_list = ['MVsolid', 'MW']
solvent_dependencies_list = []
component_dependencies_list = []
propagated_domain_notes = ''
smiles_flag = 'neutrals'

stored = {}

def calculate(solutedependencies, solventdependencies, componentdependencies, solutef, solventf, componentf):
    density = solutedependencies[0]['MW'][0] / solutedependencies[0]['MVsolid'][0]

    return round(density, round_digits), np.nan, round(0, round_digits), propagated_domain_notes, citation, units, endpoint

