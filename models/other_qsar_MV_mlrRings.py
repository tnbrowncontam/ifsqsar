import numpy
value_names = ('MVmlrRings',)
version = 1
round_digits = 2
units = 'cm^3/mol'
components = {'solute': 1, 'solvent': 0}
molecule_format = 'dev.0.0.5'
model_type = 'MLR'
intercept = False
fragmentlist = numpy.array([
                            ('[A;R1,R2]1~[A;R1,R2]~[A;R1,R2]~[A;R1,R2]~[A;R1,R2]~1', 'table ', 0.0),
                            ('[A;R1,R2]1~[A;R1,R2]~[A;R1,R2]~[A;R1,R2]~[A;R1,R2]~[A;R1,R2]~1', 'table ', 0.0),
                            ('[A;R1,R2]1~[A;R1,R2]~[A;R1,R2]~[A;R1,R2]~[A;R1,R2]~[A;R1,R2]~[A;R1,R2]~1', 'table ', 0.0),
                            ], dtype=[('smarts', 'S150'), ('description', 'S30'), ('fragstdev', float)])
smiles_flag = 'neutrals'
domain = False
coefficientarrays = numpy.array([
                                 (-3.89,),
                                 (-3.89,),
                                 (-3.89,),
                                 ], dtype=float)


def post_processing(prediction, error):
    if prediction != 0:
        prediction += 15.14
    prediction = max(-15.98, prediction)
    local_prediction = round(prediction, round_digits)
    local_error = round(error, round_digits)
    return local_prediction, local_error

