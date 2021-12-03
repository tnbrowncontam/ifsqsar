"""Meta QSAR for logSwliquid"""
import numpy as np
value_names = ('logSwliquid',)
version = 1
round_digits = 2
units = 'mg/L'
components = {'solute': 1, 'solvent': 0}
solute_dependencies_list = ['logVPliquid', 'logKaw', 'MW']
solvent_dependencies_list = []
smiles_flag = 'neutrals'


def calculate(solutedependencies, solventdependencies):
    # convert VP to WS by thermodynamic cycle
    logSw = solutedependencies['logVPliquid'][0] - np.log10(8.31446261815324) - np.log10(293.15) - solutedependencies['logKaw'][0]
    # convert from mol/m3 to mg/L
    logSw = logSw + np.log10((1000 * solutedependencies['MW'][0]) * (1/1000))
    # get average UL
    logSwUL = int(np.ceil(((solutedependencies['logVPliquid'][1]**2 + solutedependencies['logKaw'][1]**2)/2)**0.5))
    # calculate total error
    logSwerr = (solutedependencies['logVPliquid'][2]**2 + solutedependencies['logKaw'][2]**2)**0.5
    # set domain notes
    domainnotes = []
    if logSwUL <= 1:
        domainnotes.append('aggregate solute descriptor UL is in the AD')
    else:
        domainnotes.append('aggregate solute descriptor UL is out of the AD')

    return round(logSw, round_digits), logSwUL, round(logSwerr, round_digits), ', '.join(domainnotes)

