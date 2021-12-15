"""Clone of QSPR for density of liquids (densityliquid)"""
import numpy
value_names = ('densityliquid',)
version = 1
citation = 'Kotomin, A. A.; Kozlov, A. S., '\
           'Calculation of densities of organic compounds from contributions of molecular fragments. '\
           'Russ J Appl Chem 2006, 79 (6), 957-966.'
round_digits = 3
units = 'g/cm^3'
components = {'solute': 1, 'solvent': 0}
solute_dependencies_list = ['MVliquid', 'MW']
solvent_dependencies_list = []
propagated_domain_notes = ''
smiles_flag = 'neutrals'


def calculate(solutedependencies, solventdependencies):
    density = solutedependencies['MW'][0] / solutedependencies['MVliquid'][0]

    return round(density, round_digits), numpy.nan, round(0, round_digits), propagated_domain_notes, citation

