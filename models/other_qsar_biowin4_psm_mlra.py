import numpy
value_names = ('biowin4psmmlra',)
version = 1
round_digits = 2
units = 'rank'
components = {'solute': 1, 'solvent': 0}
molecule_format = 'dev.0.0.5'
model_type = 'MLRA'
intercept = False
smiles_flag = 'neutrals'
domain = False
# search for all possible configurations of 4 aromatic rings
fragmentlist = numpy.array([('c1cc2ccc3c4c2c(c1)ccc4ccc3', 'Polyaromatic hydrocarbon (4 or more rings)', 0.0),
                            ('c1ccc2-c3c4c(-c2c1)cccc4ccc3', 'Polyaromatic hydrocarbon (4 or more rings)', 0.0),
                            ('c1ccc2c(c1)cc1c(c2)ccc2c1cccc2', 'Polyaromatic hydrocarbon (4 or more rings)', 0.0),
                            ('c1ccc2c(c1)cc1c(c2)cc2c(c1)cccc2', 'Polyaromatic hydrocarbon (4 or more rings)', 0.0),
                            ('c1ccc2c(c1)ccc1c2c2ccccc2cc1', 'Polyaromatic hydrocarbon (4 or more rings)', 0.0),
                            ('c1ccc2c(c1)c1ccccc1c1c2cccc1', 'Polyaromatic hydrocarbon (4 or more rings)', 0.0),
                            ('c1ccc2c(c1)c1ccc3c(c1cc2)cccc3', 'Polyaromatic hydrocarbon (4 or more rings)', 0.0),
                            ], dtype=[('smarts', 'S150'), ('description', 'S30'), ('fragstdev', float)])
coefficientarrays = numpy.array([(1,),
                                 (1,),
                                 (1,),
                                 (1,),
                                 (1,),
                                 (1,),
                                 (1,),
                                 ], dtype=float)


def post_processing(prediction, error):
    if prediction > 0:
        local_prediction = round(-0.70224, round_digits)
    else:
        local_prediction = 0
    local_error = round(error, round_digits)
    return local_prediction, local_error

