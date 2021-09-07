"""Meta QSAR for logKow"""
import numpy
round_digits = 2
units = 'unitless'
dependencies_list = ['S', 'A', 'B', 'V', 'L']
smiles_flag = 'neutrals'


def calculate(dependencies):
    # calculate log Kow with pplfer equation from Brown 2021
    logKow = -1.57 * dependencies['S'][0] \
             - 0.16 * dependencies['A'][0] \
             - 4.05 * dependencies['B'][0] \
             + 2.71 * dependencies['V'][0] \
             + 0.41 * dependencies['L'][0] \
             + 0.38

    return round(logKow, round_digits), 0, round(0, round_digits), ''

