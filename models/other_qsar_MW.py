"""Calculates MW so that it is accessible to METAQSARs"""
import numpy
value_names = ('MW',)
version = 1
citation = 'MW'
round_digits = 6
units = 'g/mol'
components = {'solute': 1, 'solvent': 0}
molecule_format = 'v1.0.0'
model_type = 'MLR'
intercept = True
smiles_flag = 'neutrals'
domain = False
fragmentlist = numpy.array([('MW', 'molecular weight', 0.0),
                            ], dtype=[('smarts', 'S2'), ('description', 'S30'), ('fragstdev', float)])
coefficientarrays = numpy.array([(1.,),
                                 ], dtype=float)


def post_processing(prediction, error):
    local_prediction = round(prediction, round_digits)
    local_error = round(error, round_digits)
    return local_prediction, local_error

