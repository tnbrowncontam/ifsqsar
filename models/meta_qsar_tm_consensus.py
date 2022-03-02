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
components = {'solute': 1, 'solvent': 0}
solute_dependencies_list = ['tm', 'tmpplfer']
solvent_dependencies_list = []
propagated_domain_notes = ''
smiles_flag = 'neutrals'

def calculate(solutedependencies, solventdependencies):
    MP = (solutedependencies['tm'][0] + solutedependencies['tmpplfer'][0]) / 2
    MPUL = int(np.ceil(((solutedependencies['tm'][1]**2 + solutedependencies['tmpplfer'][1]**2) / 2)**0.5))
    MPerr = ((0.5 * solutedependencies['tm'][2])**2 + (0.5 * solutedependencies['tmpplfer'][2])**2)**0.5
    return round(MP, round_digits), MPUL, round(MPerr, round_digits), propagated_domain_notes, citation, units, endpoint

