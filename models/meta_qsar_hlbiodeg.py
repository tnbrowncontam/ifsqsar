"""Extrapolation of biodegradation half-life (HLbiodeg) from approximate clones of QSRPs (BIOWIN3, BIOWIN4)"""
import numpy as np
value_names = ('HLbiodeg',)
version = 1
endpoint = 'Biodegradation half-life in water (20-25degC)'
citation = 'BIOWIN3 and BIOWIN4: ' \
           'Boethling, R. S.;  Howard, P. H.;  Meylan, W.;  Stiteler, W.;  Beauman, J.; Tirado, N., '\
           'Group contribution method for predicting probability and rate of aerobic biodegradation. '\
           'Environ Sci Technol 1994, 28 (3), 459-65. '\
           'Extrapolation to HLbiodeg: '\
           'Arnot, J. A.;  Gouin, T.; Mackay, D. '\
           'Development and Application of Models of Chemical Fate in Canada - Practical methods for estimating '\
           'environmental biodegradation rates; CEMN Report No. 200503, 2005.'
round_digits = 2
units = 'hours'
chemical_inputs = {'solute min': 1, 'solute max': 1,
                   'solvent min': 0, 'solvent max': 0,
                   'component min': 0, 'component max': 0,
                   'total min': 1, 'total max': 1}
solute_dependencies_list = ['biowin3usmmlrx', 'biowin3usmmlra', 'biowin4psmmlrx', 'biowin4psmmlra']
solvent_dependencies_list = []
component_dependencies_list = []
propagated_domain_notes = ''
smiles_flag = 'neutrals'

stored = {}

def calculate(solutedependencies, solventdependencies, componentdependencies, solutef, solventf, componentf):
    domainnotes = 'generic error estimated from uncertainties of model fits'
    # estimate HLbiodeg from both biowin survey models and take the average
    biowin3 = solutedependencies[0]['biowin3usmmlrx'][0] + solutedependencies[0]['biowin3usmmlra'][0]
    biowin4 = solutedependencies[0]['biowin4psmmlrx'][0] + solutedependencies[0]['biowin4psmmlra'][0]
    HLbiodeg = ((biowin3 * -1.07 + 4.2) + (biowin4 * -1.46 + 6.51)) / 2

    return round(24*10**HLbiodeg, round_digits), np.nan, round(10**(1.11*1.96), round_digits), domainnotes, citation, units, endpoint

