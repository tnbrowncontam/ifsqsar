import numpy
value_names = ('MVliqcorr',)
version = 1
round_digits = 2
units = 'unitless'
components = {'solute': 1, 'solvent': 0}
molecule_format = 'dev.0.0.5'
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
    local_prediction = round(prediction, round_digits)
    local_error = round(error, round_digits)
    return local_prediction, local_error

