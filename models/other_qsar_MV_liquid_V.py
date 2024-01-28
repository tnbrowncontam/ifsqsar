"""QSPR for MV of liquids based on QSPR for V"""
import numpy
value_names = ('MVliquid',)
version = 2
endpoint = 'Molar volume of liquid'
citation = 'Brown T.N., Sangion A., Arnot J.A.; '\
           'Identifying Uncertainty in Physical-Chemical Property Estimation with IFSQSAR.'\
           '2024, In Prep.'
round_digits = 2
units = 'cm^3/mol'
chemical_inputs = {'solute min': 1, 'solute max': 1,
                   'solvent min': 0, 'solvent max': 0,
                   'component min': 0, 'component max': 0,
                   'total min': 1, 'total max': 1}
molecule_format = 'v1.0.0'
model_type = 'MLR'
intercept = False
smiles_flag = 'neutrals'
domain = False
fragmentlist = numpy.array([('[#5]', 'no description', 0.0),
                            ('[#6]', 'no description', 0.0),
                            ('[#7]', 'no description', 0.0),
                            ('[#8]', 'no description', 0.0),
                            ('[#9]', 'no description', 0.0),
                            ('[#14]', 'no description', 0.0),
                            ('[#15]', 'no description', 0.0),
                            ('[#16]', 'no description', 0.0),
                            ('[#17]', 'no description', 0.0),
                            ('[#32]', 'no description', 0.0),
                            ('[#33]', 'no description', 0.0),
                            ('[#34]', 'no description', 0.0),
                            ('[#35]', 'no description', 0.0),
                            ('[#53]', 'no description', 0.0),
                            ('[!#1]#[!#1]', 'no description', 0.0),
                            ('[!#1]=[!#1]', 'no description', 0.0),
                            ('[!#1]-[!#1]', 'no description', 0.0),
                            ('[!#1]:[!#1]', 'no description', 0.0),
                            ('[*H1]', 'no description', 0.0),
                            ('[*H2]', 'no description', 0.0),
                            ('[*H3]', 'no description', 0.0),
                            ('[*H4]', 'no description', 0.0),
                            ('[*H5]', 'no description', 0.0),
                            ('[*H6]', 'no description', 0.0),
                            ('sssr', 'no description', 0.0),
                            ], dtype=[('smarts', 'S11'), ('description', 'S14'), ('fragstdev', float)])
coefficientarrays = numpy.array([(33.8142,),
                                 (23.1256,),
                                 (19.2123,),
                                 (20.1725,),
                                 (24.6036,),
                                 (43.1347,),
                                 (31.6932,),
                                 (30.646,),
                                 (31.7958,),
                                 (40.8333,),
                                 (42.9797,),
                                 (34.9946,),
                                 (35.3406,),
                                 (41.606,),
                                 (-10.8407,),
                                 (-12.2985,),
                                 (-14.8916,),
                                 (-14.0663,),
                                 (4.3235,),
                                 (8.6469,),
                                 (12.9704,),
                                 (17.2938,),
                                 (21.6173,),
                                 (25.9407,),
                                 (6.708,),

                                 ], dtype=float)

stored = {'O': (round(18.01528, round_digits), 'E', numpy.nan, 'experimental value used', 'well known value', units, endpoint)}

def post_processing(prediction, error):
    local_prediction = round(prediction, round_digits)
    return local_prediction, 4.73


