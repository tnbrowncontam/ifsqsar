"""Meta QSAR for temperature of melting, consensus of two models (tmconsensus) (melting point)"""
import numpy as np
value_names = ('tmconsensus',)
version = 1
endpoint = 'Melting point - mean of QSPR and PPLFER predictions'
citation = 'tm from PPLFER: '\
           'Brown, T. N.; '\
           'Development of Iterative Fragment Selection (IFS) QSPRs for Poly-Parameter Linear Free Energy '\
           'Relationship (PPLFER) Solute Descriptors and System Parameters. '\
           'J Solution Chem 2021, In Review. '\
           'tm from QSPR: '\
           'Brown, T. N.;  Armitage, J. M.; Arnot, J. A., '\
           'Application of an Iterative Fragment Selection (IFS) Method to Estimate Entropies of Fusion and Melting'\
           'Points of Organic Chemicals. Mol Inform 2019, 38 (8-9), 1800160.'
round_digits = 2
units = 'K'
chemical_inputs = {'solute min': 1, 'solute max': 1,
                   'solvent min': 0, 'solvent max': 0,
                   'component min': 0, 'component max': 0,
                   'total min': 1, 'total max': 1}
solute_dependencies_list = ['tm', 'tmpplfer']
solvent_dependencies_list = []
component_dependencies_list = []
propagated_domain_notes = ''
smiles_flag = 'neutrals'

stored = {}

def calculate(solutedependencies, solventdependencies, componentdependencies, solutef, solventf, componentf):
    tmUL = solutedependencies[0]['tm'][1]
    if tmUL == 4:
        tmUL = 2
    tmPPUL = solutedependencies[0]['tmpplfer'][1]
    if tmPPUL == 4:
        tmPPUL = 2
    if type(tmUL) is str and not type(tmPPUL) is str:
        MP = solutedependencies[0]['tm'][0]
        MPUL = tmUL
        MPerr = solutedependencies[0]['tm'][2]
    elif not type(tmUL) is str and type(tmPPUL) is str:
        MP = solutedependencies[0]['tmpplfer'][0]
        MPUL = tmPPUL
        MPerr = solutedependencies[0]['tmpplfer'][2]
    elif type(tmUL) is str and type(tmPPUL) is str:
        MP = (solutedependencies[0]['tm'][0] + solutedependencies[0]['tmpplfer'][0]) / 2
        ulconcat = []
        if 'E' in tmUL or 'E' in tmPPUL:
            ulconcat.append('E')
        if 'U' in tmUL or 'U' in tmPPUL:
            ulconcat.append('U')
        if tmUL.replace('E', '').replace('U', '') != '' and tmPPUL.replace('E', '').replace('U', '') != '':
            ulconcat.append(str(int(np.ceil(((int(tmUL.replace('E', '').replace('U', ''))**2 + int(tmPPUL.replace('E', '').replace('U', ''))**2) / 2)**0.5))))
        elif tmUL.replace('E', '').replace('U', '') == '' and tmPPUL.replace('E', '').replace('U', '') != '':
            ulconcat.append(tmPPUL.replace('E', '').replace('U', ''))
        elif tmUL.replace('E', '').replace('U', '') != '' and tmPPUL.replace('E', '').replace('U', '') == '':
            ulconcat.append(tmUL.replace('E', '').replace('U', ''))
        MPUL = ''.join(ulconcat)
        MPerr = ((0.5 * solutedependencies[0]['tm'][2]) ** 2 + (0.5 * solutedependencies[0]['tmpplfer'][2]) ** 2) ** 0.5
    else:
        MP = (solutedependencies[0]['tm'][0] + solutedependencies[0]['tmpplfer'][0]) / 2
        if tmUL == 5 or tmPPUL == 5:
            MPUL = 5
        else:
            MPUL = int(np.ceil(((tmUL**2 + tmPPUL**2) / 2)**0.5))
        MPerr = ((0.5 * solutedependencies[0]['tm'][2])**2 + (0.5 * solutedependencies[0]['tmpplfer'][2])**2)**0.5
    return round(MP, round_digits), MPUL, round(MPerr, round_digits), propagated_domain_notes, citation, units, endpoint

