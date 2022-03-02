"""Partial clone of QSRP for MV, intended for internal use"""
import numpy
value_names = ('MVliqcorr',)
version = 1
endpoint = 'Molar volume - fraction of model for internal use'
citation = 'Kotomin, A. A.; Kozlov, A. S., '\
           'Calculation of densities of organic compounds from contributions of molecular fragments. '\
           'Russ J Appl Chem 2006, 79 (6), 957-966.'
round_digits = 2
units = 'unitless'
components = {'solute': 1, 'solvent': 0}
molecule_format = 'v1.0.0'
model_type = 'MLRA'
intercept = True
smiles_flag = 'neutrals'
domain = False
fragmentlist = numpy.array([('intercept', 'no description', 0.0),
                            ('[A]', 'aliphatic atoms present', 0.0),
                            ('[a]', 'aromatic atoms present', 0.0),
                            ], dtype=[('smarts', 'S9'), ('description', 'S30'), ('fragstdev', float)])
coefficientarrays = numpy.array([(0.90,),
                                 (-0.02,),
                                 (0.02,),
                                 ], dtype=float)


def post_processing(prediction, error):
    local_prediction = round(1/prediction, round_digits)
    local_error = round(error, round_digits)
    return local_prediction, local_error

