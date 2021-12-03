import numpy
value_names = ('MVmlrRings',)
version = 1
round_digits = 2
units = 'cm^3/mol'
components = {'solute': 1, 'solvent': 0}
molecule_format = 'dev.0.0.5'
model_type = 'MLR'
intercept = False
smiles_flag = 'neutrals'
domain = False
fragmentlist = numpy.array([
                            ('sssr', 'table ', 0.0),
                            ('C1C3CC2CC(CC1C2)C3', 'table 5', 0.0),
                            ('C12C3C4C1C5C2C3C45', 'table 5', 0.0),
                            ], dtype=[('smarts', 'S150'), ('description', 'S30'), ('fragstdev', float)])
coefficientarrays = numpy.array([
                                 (-3.89,),
                                 (3.50+(3.5-6.98),),
                                 (38.24+(38.24-33.94),)
                                 ], dtype=float)


def post_processing(prediction, error):
    if prediction != 0:
        prediction += 15.15
    prediction = max(-15.97, prediction)
    local_prediction = round(prediction, round_digits)
    local_error = round(error, round_digits)
    return local_prediction, local_error

