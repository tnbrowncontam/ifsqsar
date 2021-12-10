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
fragmentlist = numpy.array([('a1aa2aaa3a4a2a(a1)aaa4aaa3', 'Polyaromatic hydrocarbon (4 or more rings)', 0.0),
                            ('a1aaa2-a3a4a(-a2a1)aaaa4aaa3', 'Polyaromatic hydrocarbon (4 or more rings)', 0.0),
                            ('a1aaa2a(a1)aa1a(a2)aaa2a1aaaa2', 'Polyaromatic hydrocarbon (4 or more rings)', 0.0),
                            ('a1aaa2a(a1)aa1a(a2)aa2a(a1)aaaa2', 'Polyaromatic hydrocarbon (4 or more rings)', 0.0),
                            ('a1aaa2a(a1)aaa1a2a2aaaaa2aa1', 'Polyaromatic hydrocarbon (4 or more rings)', 0.0),
                            ('a1aaa2a(a1)a1aaaaa1a1a2aaaa1', 'Polyaromatic hydrocarbon (4 or more rings)', 0.0),
                            ('a1aaa2a(a1)a1aaa3a(a1aa2)aaaa3', 'Polyaromatic hydrocarbon (4 or more rings)', 0.0),
                            ('a1aa2aaaa3a2a(a1)a1aaaa2a1a3aaa2', 'Polyaromatic hydrocarbon (4 or more rings)', 0.0),
                            ('a1aaa2a(a1)aa1a(a2)aa2a1aaaa2', 'Polyaromatic hydrocarbon (4 or more rings)', 0.0),
                            ('a1aaa2a(a1)a1aa3a(a1aa2)aaaa3', 'Polyaromatic hydrocarbon (4 or more rings)', 0.0),
                            ('a1aaa2a(a1)aaa1a2a2aaaaa2a1', 'Polyaromatic hydrocarbon (4 or more rings)', 0.0),
                            ], dtype=[('smarts', 'S150'), ('description', 'S30'), ('fragstdev', float)])
coefficientarrays = numpy.array([(1,),
                                 (1,),
                                 (1,),
                                 (1,),
                                 (1,),
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

