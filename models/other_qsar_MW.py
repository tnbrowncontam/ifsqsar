"""Calculates MW so that it is accessible to METAQSARs"""
import numpy
value_names = ('MW',)
version = 1
endpoint = 'Molecular weight'
citation = 'MW'
round_digits = 6
units = 'g/mol'
chemical_inputs = {'solute min': 1, 'solute max': 1,
                   'solvent min': 0, 'solvent max': 0,
                   'component min': 0, 'component max': 0,
                   'total min': 1, 'total max': 1}
molecule_format = 'v1.0.0'
model_type = 'MLR'
intercept = True
smiles_flag = 'neutrals'
domain = False
fragmentlist = numpy.array([('MW', 'molecular weight', 0.0),
                            ], dtype=[('smarts', 'S2'), ('description', 'S30'), ('fragstdev', float)])
coefficientarrays = numpy.array([(1.,),
                                 ], dtype=float)

stored = {}

def post_processing(prediction, error):
    local_prediction = round(prediction, round_digits)
    local_error = round(error, round_digits)
    return local_prediction, local_error

