"""Meta QSAR for HLbiodeg"""
import numpy as np
value_names = ('HLbiodeg',)
version = 1
round_digits = 2
units = 'log days'
components = {'solute': 1, 'solvent': 0}
solute_dependencies_list = ['biowin3usmmlrx', 'biowin3usmmlra', 'biowin4psmmlrx', 'biowin4psmmlra']
solvent_dependencies_list = []
smiles_flag = 'neutrals'


def calculate(solutedependencies, solventdependencies):
    domainnotes = ['generic error estimated from uncertainties of model fits']
    # estimate HLbiodeg from both biowin survey models and take the average
    biowin3 = solutedependencies['biowin3usmmlrx'][0] + solutedependencies['biowin3usmmlra'][0]
    biowin4 = solutedependencies['biowin4psmmlrx'][0] + solutedependencies['biowin4psmmlra'][0]
    HLbiodeg = ((biowin3 * -1.07 + 4.2) + (biowin4 * -1.46 + 6.51)) / 2

    return round(HLbiodeg, round_digits), np.nan, round(1.11, round_digits), ', '.join(domainnotes)

