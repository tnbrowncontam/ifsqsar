import numpy
value_names = ('V',)
version = 1
round_digits = 3
units = '0.01 cm^3/mol'
components = {'solute': 1, 'solvent': 0}
molecule_format = 'old_format'
model_type = 'MLR'
intercept = False
smiles_flag = 'neutrals'
domain = False
fragmentlist = numpy.array([('[#6]', 'no description', 0.0),
                            ('[#7]', 'no description', 0.0),
                            ('[#8]', 'no description', 0.0),
                            ('[#9]', 'no description', 0.0),
                            ('[#14]', 'no description', 0.0),
                            ('[#15]', 'no description', 0.0),
                            ('[#16]', 'no description', 0.0),
                            ('[#17]', 'no description', 0.0),
                            ('[#5]', 'no description', 0.0),
                            ('[#32]', 'no description', 0.0),
                            ('[#33]', 'no description', 0.0),
                            ('[#34]', 'no description', 0.0),
                            ('[#35]', 'no description', 0.0),
                            ('[#50]', 'no description', 0.0),
                            ('[#51]', 'no description', 0.0),
                            ('[#52]', 'no description', 0.0),
                            ('[#53]', 'no description', 0.0),
                            ('[*H1]', 'no description', 0.0),
                            ('[*H2]', 'no description', 0.0),
                            ('[*H3]', 'no description', 0.0),
                            ('[*H4]', 'no description', 0.0),
                            ('[*H5]', 'no description', 0.0),
                            ('[*H6]', 'no description', 0.0),
                            ('[!#1]~[!#1]', 'no description', 0.0),
                            ], dtype=[('smarts', 'S11'), ('description', 'S14'), ('fragstdev', float)])
coefficientarrays = numpy.array([(0.1635,),
                                 (0.1439,),
                                 (0.1243,),
                                 (0.1048,),
                                 (0.2683,),
                                 (0.2487,),
                                 (0.2291,),
                                 (0.2095,),
                                 (0.1832,),
                                 (0.3102,),
                                 (0.2942,),
                                 (0.2781,),
                                 (0.2621,),
                                 (0.3935,),
                                 (0.3744,),
                                 (0.3614,),
                                 (0.3453,),
                                 (0.0215,),
                                 (0.043,),
                                 (0.0645,),
                                 (0.086,),
                                 (0.1075,),
                                 (0.129,),
                                 (-0.0656,),
                                 ], dtype=float)


def post_processing(prediction, error):
    local_prediction = round(prediction, round_digits)
    return local_prediction, error

