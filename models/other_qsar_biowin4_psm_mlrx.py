"""Partial approximate clone of QSRP for ultimate aerobic biodegradation (BIOWIN3), intended for internal use"""
import numpy
value_names = ('biowin4psmmlrx',)
version = 1
endpoint = 'BIOWIN 4 - fraction of model for internal use'
citation = 'Boethling, R. S.;  Howard, P. H.;  Meylan, W.;  Stiteler, W.;  Beauman, J.; Tirado, N., '\
           'Group contribution method for predicting probability and rate of aerobic biodegradation. '\
           'Environ Sci Technol 1994, 28 (3), 459-65.'
round_digits = 2
units = 'rank'
components = {'solute': 1, 'solvent': 0}
molecule_format = 'v1.0.0'
model_type = 'MLRX'
intercept = True
smiles_flag = 'neutrals'
domain = False
fragmentlist = numpy.array([('[Nv3;H1&$(N-[#6!$(C=O)]),H0&$(N(-[#6!$(C=O)])-[#6!$(C=O)])]-[NH0v3]=O', 'Nitroso [-N-N=O]', 0.0),
                            ('[C!$(C=O);D2,D1]~[CD2]~[CD2]-[CH3]', 'Linear C4 terminal chain [CCC-CH3]', 0.0),
                            ('[OH$(O-[C!$(C=[!#6])])!$(O-[CH0X4])]', 'Aliphatic alcohol [-OH]', 0.0),
                            ('[OH$(O-c)]', 'Aromatic alcohol [-OH]', 0.0),
                            ('[OH]-[C$(C-C)]=O', 'Aliphatic acid [-C(=O)-OH]', 0.0),
                            ('[OH]-[C$(C-c)]=O', 'Aromatic acid [-C(=O)-OH]', 0.0),
                            ('O=[C!$(C-O);H1,H2]', 'Aldehyde [-CHO]', 0.0),
                            ('O=[C!$(C(O)[O,N,S])]-[O$([#8]-,:[#6!$(C=O)])]', 'Ester [-C(=O)-O-C]', 0.0),
                            ('[Nv3;H2,H1&$(N-[P,#6;!$(C=O)]),H0&$(N(-[#6!$(C=O)])-[#6!$(C=O)]),H0&$(N=[#6])]-[C!$(C(-N)-[O,N,S])!$(C(=N)-[O,N,S])]=[S,O]', 'Amide [-C(=O)-N or -C(=S]-N]', 0.0),
                            ('[Nv3;H1,H0&$(N-[P,#6;!$(C=O)])]=[C!$(C(=N)(-[H1;S,O])-[O,N,S])]-[H1;S,O]', 'Amide [-C(=O)-N or -C(=S]-N]', 0.0),
                            ('c1ncncn1', 'Triazine ring (symmetric)', 0.0),
                            ('[Cl$(Cl-C)]', 'Aliphatic chloride [-CL]', 0.0),
                            ('[Cl$(Cl-c)]', 'Aromatic chloride [-CL]', 0.0),
                            ('[Br$(Br-C)]', 'Aliphatic bromide [-Br]', 0.0),
                            ('[Br$(Br-c)]', 'Aromatic bromide [-Br]', 0.0),
                            ('[I$(I-c)]', 'Aromatic iodide [-I]', 0.0),
                            ('[F$(F-c)]', 'Aromatic fluoride [-F]', 0.0),
                            ('[CH0X4!$(C(-F)(-F)(-F)-[!#1!#9])]', 'Carbon with 4 single bonds & nohydrogens', 0.0),
                            ('[N$(N-c)](=O)=O', 'Aromatic nitro [-NO2]', 0.0),
                            ('[Nv3;H2&$(N-C),H1&$(N(-C)-C)]', 'Aliphatic amine [-NH2  or  -NH-]', 0.0),
                            ('[Nv3$(N-c)!$(N-C=[O,S,N]);H2,H1]', 'Aromatic amine [-NH2  or  -NH-]', 0.0),
                            ('[Nv3$(N=c)!$(N-C=[O,S,N]);H1,H0]', 'Aromatic amine [-NH2  or  -NH-]', 0.0),
                            ('[C$(C-[#6])]#[Nv3]', 'Cyanide / Nitriles [-C#N]', 0.0),
                            ('[S$(S-c)](=O)(=O)[OH]', 'Sulfonic acid / salt -> aromatic attach', 0.0),
                            ('[S$(S-C)](=O)(=O)[OH]', 'Sulfonic acid / salt -> aliphatic attach', 0.0),
                            ('n1[cR1][cR1][cR1][cR1][cR1]1', 'Pyridine ring', 0.0),
                            ('[O$(O(-c)-[#6!$(C=O)])]', 'Aromatic ether [-O-aromatic carbon]', 0.0),
                            ('[O$(O(-C)-C)!$(O-C=O)]', 'Aliphatic ether [C-O-C]', 0.0),
                            ('[C$(C(-[#6!$([#6]=O)])-[#6!$([#6]=O)])]=O', 'Ketone [-C-C(=O)-C-]', 0.0),
                            ('[Nv3$(N(-[#6])(-[#6])-[#6])!$(N-C=[O,S,N])]', 'Tertiary amine', 0.0),
                            ('[O,S]=P(-[O,S;$(*-[#6])])(-[O,S;$(*-[#6])])-[O,S;$(*-[#6])]', 'Phosphate ester', 0.0),
                            ('[C!H0$(C-c)!$(C=,#*)!$(C(F)(F)F)!$(C-[O,S]-[C,P]=[O,S])!$(C-C#N)!$(C-[OH])]', 'Alkyl substituent on aromatic ring', 0.0),
                            ('[Nv3$(N-[#6])]=[Nv3$(N-[#6])]', 'Azo group [-N=N-]', 0.0),
                            ('[O,S;v2&$(*(-C)-[#6])]-C(=[O,S])-[Nv3!$(N~[!#1!#6])]', 'Carbamate or Thiocarbamate', 0.0),
                            ('[C$(C(-F)(-F)(-F)-[!#1!#9])]', 'Trifluoromethyl group [-CF3]', 0.0),
                            ('[cH0$(c-[!#1])]1[cH1][cH1][cH1][cH1][cH1]1', 'Unsubstituted phenyl group (C6H5-)', 0.0),
                            ('[cH]1:[cH]:[cH]:[cH]:[cH]:[cH]:1', 'Unsubstituted aromatic [3 or less rings): benzene', 0.0),
                            ('[cH]1:[cH]:[cH]:[cH]:[cH0R2x3]2:[cH0R2x3]:1:[cH]:[cH]:[cH]:[cH]:2', 'Unsubstituted aromatic [3 or less rings): naphthalene', 0.0),
                            ('[cH]1:[cH]:[cH]:[cH0R2x3]2-[cH0R2x3](:[cH]:[cH]:1):[cH]:[cH]:[cH]:2', 'Unsubstituted aromatic [3 or less rings): 7:5 naphthalene', 0.0),
                            ('[cH]2:[cH]:[cH]:[cH]:[cH0R2x3]3:[cH]:[cH0R2x3]1:[cH]:[cH]:[cH]:[cH]:[cH0R2x3]:1:[cH]:[cH0R2x3]:2:3', 'Unsubstituted aromatic [3 or less rings): anthracene', 0.0),
                            ('[cH]2:[cH]:[cH]:[cH]:[cH0R2x3]3:[cH]:[cH]:[cH0R2x3]1:[cH]:[cH]:[cH]:[cH]:[cH0R2x3]:1:[cH0R2x3]:2:3', 'Unsubstituted aromatic [3 or less rings): phenanthrene', 0.0),
                            ('MW', 'molecular weight', 0.0),
                            ('intercept', '', 0.0),
                            ], dtype=[('smarts', 'S150'), ('description', 'S30'), ('fragstdev', float)])
coefficientarrays = numpy.array([(0.01848,),
                                 (0.26907,),
                                 (0.12945,),
                                 (0.03969,),
                                 (0.38557,),
                                 (0.00775,),
                                 (0.19664,),
                                 (0.22896,),
                                 (0.20543,),
                                 (0.20543,),
                                 (-0.05752,),
                                 (-0.10061,),
                                 (-0.16534,),
                                 (0.03538,),
                                 (-0.15351,),
                                 (-0.12707,),
                                 (0.01346,),
                                 (-0.15344,),
                                 (-0.10838,),
                                 (0.04328,),
                                 (-0.10838,),
                                 (-0.10838,),
                                 (-0.0652,),
                                 (0.02162,),
                                 (0.17714,),
                                 (-0.01874,),
                                 (0.07712,),
                                 (-0.00974,),
                                 (-0.02222,),
                                 (-0.288,),
                                 (0.46535,),
                                 (-0.06853,),
                                 (-0.05279,),
                                 (0.19363,),
                                 (-0.2744,),
                                 (0.00489,),
                                 (-0.3428,),
                                 (-0.3428,),
                                 (-0.3428,),
                                 (-0.3428,),
                                 (-0.3428,),
                                 (-0.001442756,),
                                 (3.847737,),
                                 ], dtype=float)


def post_processing(prediction, error):
    local_prediction = round(prediction, round_digits)
    local_error = round(error, round_digits)
    return local_prediction, local_error

