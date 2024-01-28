"""QSPR for MV of solids based on QSPR for V"""
import numpy
value_names = ('MVsolid',)
version = 2
endpoint = 'Molar volume of solid'
citation = 'Brown T.N., Sangion A., Arnot J.A.; '\
           'Identifying Uncertainty in Physical-Chemical Property Estimation with IFSQSAR.'\
           '2024, In Prep.'
round_digits = 2
units = 'cm^3/mol'
chemical_inputs = {'solute min': 1, 'solute max': 1,
                   'solvent min': 0, 'solvent max': 0,
                   'component min': 0, 'component max': 0,
                   'total min': 1, 'total max': 1}
solute_dependencies_list = ['MVliquid']
solvent_dependencies_list = []
component_dependencies_list = []
propagated_domain_notes = ''
smiles_flag = 'neutrals'

stored = {}

def calculate(solutedependencies, solventdependencies, componentdependencies, solutef, solventf, componentf):
    # Apply solid -> liquid correction defined in Kotomin and Kozolov 2006
    MV = solutedependencies[0]['MVliquid'][0] * 0.962

    return round(MV, round_digits), numpy.nan, 12.64, propagated_domain_notes, citation, units, endpoint

