"""Extrapolation of biodegradation half-life (HLbiodeg) from approximate clones of QSRPs (BIOWIN3, BIOWIN4)"""
import numpy as np
value_names = ('HLbiodeg',)
version = 1
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
components = {'solute': 1, 'solvent': 0}
solute_dependencies_list = ['biowin3usmmlrx', 'biowin3usmmlra', 'biowin4psmmlrx', 'biowin4psmmlra']
solvent_dependencies_list = []
propagated_domain_notes = ''
smiles_flag = 'neutrals'


def calculate(solutedependencies, solventdependencies):
    domainnotes = 'generic error estimated from uncertainties of model fits'
    # estimate HLbiodeg from both biowin survey models and take the average
    biowin3 = solutedependencies['biowin3usmmlrx'][0] + solutedependencies['biowin3usmmlra'][0]
    biowin4 = solutedependencies['biowin4psmmlrx'][0] + solutedependencies['biowin4psmmlra'][0]
    HLbiodeg = ((biowin3 * -1.07 + 4.2) + (biowin4 * -1.46 + 6.51)) / 2

    return round(24*10**HLbiodeg, round_digits), np.nan, round(10**(1.11*1.96), round_digits), domainnotes, citation, units

