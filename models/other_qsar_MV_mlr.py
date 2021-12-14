"""Partial clone of QSRP for MV, intended for internal use"""
import numpy
value_names = ('MVmlr',)
version = 1
citation = 'Kotomin, A. A.; Kozlov, A. S., '\
           'Calculation of densities of organic compounds from contributions of molecular fragments. '\
           'Russ J Appl Chem 2006, 79 (6), 957-966.'
round_digits = 2
units = 'cm^3/mol'
components = {'solute': 1, 'solvent': 0}
molecule_format = 'dev.0.0.5'
model_type = 'MLR'
intercept = False
smiles_flag = 'neutrals'
domain = False
fragmentlist = numpy.array([
                            ('[OH1D1]-[CR0]~[CR0]-[OH1D1]', 'table 2', 0.0),
                            ('[CR0](-[NH0D3R0v5](=O)=O)~[CR0]-[NH0D3R0v5](=O)=O', 'table 2', 0.0),
                            ('[CR0](-[OH0D2R0]-[NH0D3R0v5](=O)=O)~[CR0]-[NH0D3R0v5](=O)=O', 'table 2', 0.0),
                            ('[CR0](-[OH0D2R0]-[NH0D3R0v5](=O)=O)~[CR0]-[OH0D2R0]-[NH0D3R0v5](=O)=O', 'table 2', 0.0),
                            ('[CR0](-[NH0D3R0v3]-[NH0D3R0v5](=O)=O)~[CR0]-[CH0D3R0]=O', 'table 2', 0.0),
                            ('[CR0](-[NH0D3R0v5](=O)=O)~[CR0]-[CH0D3R0]=O', 'table 2', 0.0),
                            ('[CR0](-[NH1D2R0v3]-[NH0D3R0v5](=O)=O)-[NH1D2R0v3]-[NH0D3R0v5](=O)=O', 'table 2', 0.0),
                            ('[CR0](-[NH1D2R0v3]-[NH0D3R0v5](=O)=O)-[NH0D3R0v3]-[NH0D3R0v5](=O)=O', 'table 2', 0.0),
                            ('[CR0](-[NH0D3R0v3]-[NH0D3R0v5](=O)=O)-[NH0D3R0v3]-[NH0D3R0v5](=O)=O', 'table 2', 0.0),
                            ('[CR0](=[NH0D2R0v3]-[NH0D3R0v5](=O)=O)-[NH1D2R0v3]-[NH0D3R0v5](=O)=O', 'table 2', 0.0),
                            ('[CR0](-[NH0D3R0v3]-[NH0D3R0v5](=O)=O)-[NH0D2R0v3]=[NH0D2R0v5]#[NH0D1R0v3]', 'table 2', 0.0),
                            ('[CR0](=[NH0D2R0v3]-[NH0D3R0v5](=O)=O)-[NH0D3R0v3]-[NH0D3R0v5](=O)=O', 'table 2', 0.0),
                            ('[CR0](-[NH0D3R0v3](-F)-F)-[NH0D3R0v5](=O)=O', 'table 2', 0.0),
                            ('[CR0](-Cl)-Cl', 'table 2', 0.0),
                            ('[CR0](-Cl)-[NH0D3R0v5](=O)=O', 'table 2', 0.0),
                            ('[CR0](-[OH1D1])~[CR0]-[CH0D3R0](=O)-[OH1]', 'table 2', 0.0),
                            ('[CR0](-[CH0D3R0]=O)~[CR0]-[CH0D3R0]=O', 'table 2', 0.0),
                            ('[CR1](-[NH0D3R0v5$(*(-[CR])-[CR])](=O)=O)-@[CR2]-@-[CR1]-[NH0D3R0v5](=O)=O', 'table 4', 0.0),
                            ('[NR](-[CH0D3R$(*(-[NR])-[NR])]=O)~[C,N;R]-[NH1D2Rv3$(*(-[C,N;R])-[NR])]', 'table 4', 0.0),
                            ('[NR](-[CH0D3R$(*(-[NR])-[NR])]=O)~[CR]-[NH0D3Rv3$(*(-[CR])-[CR])]-[NH0D3R0v5](=O)=O', 'table 4', 0.0),
                            ('[cH1]', 'table 6', 0.0),
                            ('[nH1]', 'table 6', 0.0),
                            ('[c$(c1:c:c:c:c:c:1)]-!@[c$(c1:c:c:c:c:c:1)]', 'table 7', 0.0),
                            ('c~!@[c!$(c1:c:c:c:c:c:1)]', 'table 7', 0.0),
                            ('[c$(c1:c:c:c:c:c:1)]-!@n', 'table 7', 0.0),
                            ('[c!$(c1:c:c:c:c:c:1)]-!@n', 'table 7', 0.0),
                            ], dtype=[('smarts', 'S150'), ('description', 'S30'), ('fragstdev', float)])
coefficientarrays = numpy.array([
                                 (1.56,),
                                 (0.42,),
                                 (0.42,),
                                 (0.42,),
                                 (3.88,),
                                 (3.88,),
                                 (1.24,),
                                 (1.24,),
                                 (1.24,),
                                 (1.24,),
                                 (1.78,),
                                 (11.06,),
                                 (1.38,),
                                 (2.23,),
                                 (2.23,),
                                 (4.42,),
                                 (4.50,),
                                 (8.18,),
                                 (-2.35,),
                                 (2.29,),
                                 (7.96,),
                                 (5.63,),
                                 (5.39,),
                                 (1.97,),
                                 (2.28,),
                                 (1.16,),
                                 ], dtype=float)


def post_processing(prediction, error):
    local_prediction = round(prediction, round_digits)
    local_error = round(error, round_digits)
    return local_prediction, local_error

