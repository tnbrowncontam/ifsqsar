"""Meta QSAR for logKsa"""
import numpy as np
value_names = ('logKsa',)
version = 1
endpoint = 'Log of solvent-air partition coefficient - user-defined solvent'
citation = 'Brown, T. N.; '\
           'QSPRs for Predicting Equilibrium Partitioning in Solvent-Air Systems '\
           'from the Chemical Structures of Solutes and Solvents. '\
           'J Solution Chem 2022, (https://doi.org/10.1007/s10953-022-01162-2).'
round_digits = 2
units = 'log L[a]/L[s]'
chemical_inputs = {'solute min': 1, 'solute max': 1,
                   'solvent min': 1, 'solvent max': 1,
                   'component min': 0, 'component max': 0,
                   'total min': 2, 'total max': 2}
solute_dependencies_list = ['E', 'S', 'A', 'B', 'V', 'L']
solvent_dependencies_list = ['E', 'S', 'A', 'B', 'V', 'L', 's', 'a', 'b', 'v', 'l', 'c', 'state']
component_dependencies_list = []
propagated_domain_notes = ''
smiles_flag = 'neutrals'

emptrainset = np.array([[0, 0, 0, 0, 2.363, 7.714],
                        [0, 0, 0, 0, 1.6585, 5.191],
                        [0, 0, 0, 0, 1.2358, 3.677],
                        [0, 0, 0, 0, 1.5176, 4.686],
                        [0, 0, 0, 0, 1.2358, 3.106],
                        [0, 0.25, 0, 0.45, 1.2945, 3.924],
                        [-0.06, 0.17, 0, 0.57, 1.0127, 2.501],
                        [0.024, 0.22, 0, 0.59, 0.8718, 2.38],
                        [0.041, 0.25, 0, 0.45, 0.7309, 2.015],
                        [0.07, 0.6, 0, 0.45, 1.0284, 3.353],
                        [0.21, 0.4, 0, 0.1, 0.7946, 2.722],
                        [0.09, 0.6, 0, 0.45, 0.8875, 2.819],
                        [0.66, 0.56, 0, 0.16, 0.9982, 3.939],
                        [0.61, 0.51, 0, 0.15, 0.9982, 3.778],
                        [0.61, 0.52, 0, 0.16, 0.9982, 3.839],
                        [0.46, 0.38, 0, 0, 0.7391, 2.823],
                        [0.23, 1.09, 0, 0.76, 1.4922, 5.618],
                        [0.62, 0.52, 0, 0.16, 0.9982, 3.839],
                        [0.191, 0.42, 0.37, 0.48, 1.5763, 5.61],
                        [0.06, 0.58, 0, 0.53, 0.9462, 3.412],
                        [0.14, 0.64, 0, 0.45, 0.6057, 1.911],
                        [0.2, 0.42, 0.37, 0.48, 1.2945, 4.619],
                        [0.21, 0.42, 0.37, 0.48, 1.0127, 3.61],
                        [0.17, 0.33, 0.33, 0.56, 1.0127, 3.179],
                        [0.88, 0.73, 0, 0.09, 0.8914, 4.041],
                        [0.21, 0.42, 0.37, 0.48, 1.1536, 4.115],
                        [0.17, 0.7, 0, 0.51, 0.6879, 2.287],
                        [0.72, 0.65, 0, 0.07, 0.8388, 3.657],
                        [0.82, 1.01, 0, 0.48, 1.0139, 4.501],
                        [0.22, 0.42, 0.37, 0.48, 0.8718, 3.106],
                        [0.19, 0.9, 0, 0.36, 0.686, 2.548],
                        [0.74, 1.11, 0, 0.33, 0.8711, 4.039],
                        [0.29, 0.52, 0, 0.48, 0.6223, 2.636],
                        [0.142, 0.54, 0, 0.57, 0.6644, 2.328],
                        [0.201, 0.53, 0.26, 0.83, 1.0714, 3.656],
                        [0.22, 0.39, 0.37, 0.48, 0.7309, 2.413],
                        [0.4, 0.86, 0, 0.56, 0.8611, 3.792],
                        [0.2, 0.48, 0.21, 0.91, 0.9305, 3.214],
                        [0.19, 0.39, 0.37, 0.48, 0.8718, 3.011],
                        [0.676, 0.43, 0, 0.12, 0.5077, 2.106],
                        [0.21, 0.5, 0.3, 0.83, 0.9305, 3.31],
                        [0.22, 0.42, 0.37, 0.48, 0.7309, 2.601],
                        [0.33, 0.75, 0, 0.64, 0.681, 2.892],
                        [0.363, 1.38, 0, 0.8, 0.7877, 3.639],
                        [0.22, 0.36, 0.33, 0.56, 0.7309, 2.338],
                        [0.42, 0.64, 0.1, 0.11, 0.6352, 2.573],
                        [0.237, 0.55, 0.29, 0.82, 0.7896, 2.719],
                        [0.18, 0.7, 0.04, 0.49, 0.547, 1.696],
                        [0.53, 0.93, 0, 0.84, 0.9609, 3.971],
                        [0.49, 1.3, 0, 0.79, 0.82, 3.832],
                        [0.56, 1.2, 0, 0.92, 0.8445, 3.905],
                        [0.21, 0.36, 0.33, 0.56, 0.59, 1.764],
                        [0.56, 1.38, 0, 0.99, 0.8787, 4.25],
                        [0.25, 0.42, 0.37, 0.48, 0.4491, 1.485],
                        [0.24, 0.42, 0.37, 0.48, 0.59, 2.031],
                        [0.16, 0.9, 0.02, 0.36, 0.5451, 2.082],
                        [0.367, 1.31, 0, 0.74, 0.6468, 3.173],
                        [0.96, 0.96, 0.26, 0.41, 0.8162, 3.934],
                        [0, 0.6, 0.59, 0.46, 0.1673, 0.245],
                        [0.803, 0.87, 0.39, 0.56, 0.916, 4.221],
                        [0.38, 1.33, 0.54, 1.42, 1.1888, 5.13],
                        [0.34, 1.3, 0.37, 0.71, 0.7877, 3.844],
                        [0.39, 0.88, 0.52, 1.15, 0.8483, 3.57],
                        [0.35, 1.28, 0.4, 0.71, 0.6468, 2.974],
                        [0.47, 1.31, 0.64, 0.57, 0.365, 2.447],
                        [0.4, 0.9, 0.58, 0.78, 0.5078, 2.661]
                        ], dtype=float)
xtxi = np.linalg.inv(np.matmul(emptrainset.T, emptrainset))

stored = {}


def aggregateUL(ULlist, roundup=False):
    ise = False
    isu = False
    is5 = False
    numul = None
    numcount = 0
    for UL in ULlist:
        if type(UL) is str:
            if UL == 'E':
                ise = True
            elif UL == 'U':
                isu = True
            elif UL == 'EU':
                ise = True
                isu = True
            else:
                if 'E' in UL:
                    ise = True
                if 'U' in UL:
                    isu = True
                if UL.replace('E', '').replace('U', '') == '5':
                    is5 = True
                if numul is None:
                    numul = float(UL.replace('E', '').replace('U', ''))**2
                else:
                    numul += float(UL.replace('E', '').replace('U', ''))**2
                numcount += 1
        else:
            if UL == 5:
                is5 = True
            if numul is None:
                numul = UL**2
            else:
                numul += UL**2
            numcount += 1
    ulconcat = []
    if ise:
        ulconcat.append('E')
    if isu:
        ulconcat.append('U')
    if numul is not None:
        if is5:
            ULtot = 5
        elif roundup:
            ULtot = int(np.ceil((numul/numcount)**0.5))
        else:
            ULtot = (numul/numcount)**0.5
    else:
        ULtot = None
    if len(ulconcat) and ULtot is not None:
        ulconcat.append(str(ULtot))
        ULtot = ''.join(ulconcat)
    elif len(ulconcat) and ULtot is None:
        ULtot = ''.join(ulconcat)
    return ULtot


def calculate(solutedependencies, solventdependencies, componentdependencies, solutef, solventf, componentf):
    # determine if the solvent is a liquid
    domainnotes = [propagated_domain_notes]
    # calculate leverage of the solvent vs. the empirical correlations training dataset
    x = np.array([solventdependencies[0]['E'][0],
                  solventdependencies[0]['S'][0],
                  solventdependencies[0]['A'][0],
                  solventdependencies[0]['B'][0],
                  solventdependencies[0]['V'][0],
                  solventdependencies[0]['L'][0]])
    leverage = np.matmul(np.matmul(x, xtxi), x.T)
    if leverage < 1.5 * 6 / 66:
        ULemptra = 0
    elif leverage < 3 * 6 / 66:
        ULemptra = 1
    elif leverage > 1:
        ULemptra = 3
    else:
        ULemptra = 2
    # get state of solvent and set error scaling
    if solventdependencies[0]['state'][3] == 'likely liquid':
        domainnotes.append('solvent phase: likely liquid, prediction error scaling = 1')
        phaseerrorscaling = 1
    elif solventdependencies[0]['state'][3] == 'maybe liquid':
        domainnotes.append('solvent phase: maybe liquid, prediction error scaling = 1.25')
        phaseerrorscaling = 1.25
    elif solventdependencies[0]['state'][3] == 'likely solid':
        domainnotes.append('solvent phase: likely solid, prediction error scaling = 2')
        phaseerrorscaling = 2
    elif solventdependencies[0]['state'][3] == 'maybe solid':
        domainnotes.append('solvent phase: maybe solid, prediction error scaling = 1.5')
        phaseerrorscaling = 1.5
    elif solventdependencies[0]['state'][3] == 'likely gas':
        domainnotes.append('solvent phase: likely gas, prediction error scaling = 1.5')
        phaseerrorscaling = 1.5
    elif solventdependencies[0]['state'][3] == 'maybe gas':
        domainnotes.append('solvent phase: maybe gas, prediction error scaling = 1.25')
        phaseerrorscaling = 1.25
    else:
        domainnotes.append('solvent phase: unclassified, prediction error scaling = 2')
        phaseerrorscaling = 2
    # calculate aggregate UL of solute descriptors
    ULlist = []
    for sltdes in ['S', 'A', 'B', 'L']:
        if solutedependencies[0][sltdes][1] in ['E', 'U']:
            ULlist.append(solutedependencies[0][sltdes][1])
        elif solutedependencies[0][sltdes][1] < 4:
            ULlist.append(solutedependencies[0][sltdes][1])
        elif solutedependencies[0][sltdes][1] == 4 and sltdes != 'L':
            ULlist.append(1)
        elif solutedependencies[0][sltdes][1] == 4 and sltdes == 'L':
            ULlist.append(2)
        elif solutedependencies[0][sltdes][1] == 5:
            ULlist.append(5)
        elif solutedependencies[0][sltdes][1] == 6:
            ULlist.append(3)
    ULslt = aggregateUL(ULlist)
    # calculate aggregate UL of system parameters from QSPRs
    ULlist = []
    for slvpar in ['s', 'a', 'b', 'v', 'l', 'c']:
        if solventdependencies[0][slvpar][1] in ['E', 'U']:
            ULlist.append(solventdependencies[0][slvpar][1])
        elif solventdependencies[0][slvpar][1] < 4:
            ULlist.append(solventdependencies[0][slvpar][1])
        elif solventdependencies[0][slvpar][1] == 4 and slvpar in ['s', 'a', 'b']:
            ULlist.append(1)
        elif solventdependencies[0][slvpar][1] == 4 and slvpar in ['v', 'l', 'c']:
            ULlist.append(2)
        elif solventdependencies[0][slvpar][1] == 5:
            ULlist.append(5)
        elif solventdependencies[0][slvpar][1] == 6:
            ULlist.append(3)
    ULslvqspr = aggregateUL(ULlist)
    # calculate relative solute descriptors, UL and error
    # Vri
    Vri = 1 / solventdependencies[0]['V'][0]
    Vul = 0
    Ver = 0
    # Eri
    Eri = Vri * solventdependencies[0]['E'][0]
    Eul = solventdependencies[0]['E'][1]
    if Eul == 4:
        Eul = 1
    elif Eul == 6:
        Eul = 3
    Eer = Vri * solventdependencies[0]['E'][2]
    # Sri
    Sri = Vri * solventdependencies[0]['S'][0]
    Sul = solventdependencies[0]['S'][1]
    if Sul == 4:
        Sul = 1
    elif Sul == 6:
        Sul = 3
    Ser = Vri * solventdependencies[0]['S'][2]
    # Ari
    Ari = Vri * solventdependencies[0]['A'][0]
    Aul = solventdependencies[0]['A'][1]
    if Aul == 4:
        Aul = 1
    elif Aul == 6:
        Aul = 3
    Aer = Vri * solventdependencies[0]['A'][2]
    # Bri
    Bri = Vri * solventdependencies[0]['B'][0]
    Bul = solventdependencies[0]['B'][1]
    if Bul == 4:
        Bul = 1
    elif Bul == 6:
        Bul = 3
    Ber = Vri * solventdependencies[0]['B'][2]
    # Lri
    Lri = Vri * solventdependencies[0]['L'][0]
    Lul = solventdependencies[0]['L'][1]
    if Lul == 4:
        Lul = 1
    elif Lul == 6:
        Lul = 3
    Ler = Vri * solventdependencies[0]['L'][2]
    # IAB
    IAB = (Ari * Bri)**0.5
    IABul = aggregateUL([Aul, Bul])
    IABer = 0
    if IAB > 0:
        IABer = Ari * Bri * ((Aer/Ari)**2 + (Ber/Bri)**2)**0.5
        IABer = IAB * 0.5 * IABer / (Ari * Bri)
    # IVL
    IVL = -np.log((np.exp(-Vri) + np.exp(-Lri)) / 2)
    IVLul = Lul
    IVLer = 0.5 * abs(np.exp(-Lri)) * Ler / (0.5 * (np.exp(-Vri) + np.exp(-Lri)))
    # calculate empirical system parameters, UL and error
    # s emp
    if Eri <= 0 and Sri <= 0:
        semp = 0
        if type(Eul) is str or type(Sul) is str:
            sempul = aggregateUL([Eul, Sul])
        elif Eul == 5 or Sul == 5:
            sempul = 5
        else:
            sempul = 1
        semperr = 0
    elif Eri <= 0.15 and Sri <= 0.15:
        semp = 0.16
        if type(Eul) is str or type(Sul) is str:
            sempul = aggregateUL([Eul, Sul])
        elif Eul == 5 or Sul == 5:
            sempul = 5
        else:
            sempul = 1
        semperr = semp / 1.96
    else:
        semp = 1.327 * Sri - 0.405 * Vri - 0.239 * Lri - 0.507 * IAB + 1.871
        sempul = aggregateUL([Sul, Lul, IABul])
        if semp < 0:
            semp = 0
            if sempul not in [5, 'E5', 'U5', 'EU5']:
                sempul = 3
        semperr = -0.405**2 * Vri**2 * 0.010**2 + (-0.239 * Lri)**2 * ((0.007 / -0.239)**2 + (Ler / Lri)**2) + 0.028**2
        if Sri != 0:
            semperr += (1.327 * Sri)**2 * ((0.011 / 1.327)**2 + (Ser / Sri)**2)
        if IAB != 0:
            semperr += (-0.507 * IAB)**2 * ((0.008 / -0.507)**2 + (IABer / IAB)**2)
        semperr = semperr ** 0.5
    # a emp
    if Bri == 0:
        aemp = 0
        if type(Bul) is str:
            aempul = Bul
        elif Bul == 5:
            aempul = 5
        else:
            aempul = 1
        aemperr = 0
    elif Bri < 0.1:
        aemp = 0.54
        if type(Bul) is str:
            aempul = Bul
        elif Bul == 5:
            aempul = 5
        else:
            aempul = 1
        aemperr = aemp / 1.96
    else:
        aemp = 2.845 * Bri + 0.421 * IAB - 1.582 * Sri + 0.607 * Lri
        aempul = aggregateUL([Bul, IABul, Sul, Lul])
        if aemp < 0:
            aemp = 0
            if aempul not in [5, 'E5', 'U5', 'EU5']:
                aempul = 3
        aemperr = (0.607 * Lri)**2 * ((0.004 / 0.607)**2 + (Ler / Lri)**2)
        if Bri != 0:
            aemperr += (2.845 * Bri)**2 * ((0.030 / 2.845)**2 + (Ber / Bri)**2)
        if IAB != 0:
            aemperr += (0.421 * IAB)**2 * ((0.021 / 0.421)**2 + (IABer / IAB)**2)
        if Sri != 0:
            aemperr += (-1.582 * Sri)**2 * ((0.016 / -1.582)**2 + (Ser / Sri)**2)
        aemperr = aemperr ** 0.5
    # b emp
    if Ari == 0:
        bemp = 0
        if type(Aul) is str:
            bempul = Aul
        elif Aul == 5:
            bempul = 5
        else:
            bempul = 1
        bemperr = 0
    else:
        bemp = 0.747 * Ari + 0.378 * IAB - 0.228 * Sri + 0.312 * Vri
        bempul = aggregateUL([Aul, IABul, Sul])
        bemperr = 0.312**2 * Vri**2 * 0.007**2
        if bemp < 0:
            bemp = 0
            if bempul not in [5, 'E5', 'U5', 'EU5']:
                bempul = 3
        if Ari != 0:
            bemperr += (0.747 * Ari)**2 * ((0.018 / 0.747)**2 + (Aer / Ari)**2)
        if IAB != 0:
            bemperr += (0.378 * IAB)**2 * ((0.026 / 0.378)**2 + (IABer / IAB)**2)
        if Sri != 0:
            bemperr += (-0.228 * Sri)**2 * ((0.006 / -0.228)**2 + (Ser / Sri)**2)
        bemperr = bemperr ** 0.5
    # v emp
    vemp = -0.583 * Eri - 0.279 * Sri + 0.980 * IVL - 0.918 * IAB - 0.451
    vempul = aggregateUL([Eul, Sul, IVLul, IABul])
    vemperr = (0.980 * IVL)**2 * ((0.012 / 0.980)**2 + (IVLer / IVL)**2) + 0.018**2
    if Eri != 0:
        vemperr += (-0.583 * Eri)**2 * ((0.016 / -0.583)**2 + (Eer / Eri)**2)
    if Sri != 0:
        vemperr += (-0.279 * Sri)**2 * ((0.007 / -0.279)**2 + (Ser / Sri)**2)
    if IAB != 0:
        vemperr += (-0.918 * IAB)**2 * ((0.008 / -0.918)**2 + (IABer / IAB)**2)
    vemperr = vemperr**0.5
    # l emp
    lemp = 0.215 * Eri - 0.119 * Sri - 0.152 * IVL + 0.053 * IAB + 0.986
    lempul = aggregateUL([Eul, Sul, IVLul, IABul])
    lemperr = (-0.152 * IVL)**2 * ((0.003 / -0.152)**2 + (IVLer / IVL)**2) + 0.004**2
    if Eri != 0:
        lemperr += (0.215 * Eri)**2 * ((0.004 / 0.215)**2 + (Eer / Eri)**2)
    if Sri != 0:
        lemperr += (-0.119 * Sri)**2 * ((0.002 / -0.119)**2 + (Ser / Sri)**2)
    if IAB != 0:
        lemperr += (0.053 * IAB)**2 * ((0.002 / 0.053)**2 + (IABer / IAB)**2)
    lemperr = lemperr**0.5
    # c emp
    cemp = -0.253 * Eri - 0.296 * Bri + 0.157 * IVL - 0.189
    cempul = aggregateUL([Eul, Bul, IVLul])
    cemperr = (0.157 * IVL)**2 * ((0.004 / 0.157)**2 + (IVLer / IVL)**2) + 0.005**2
    if Eri != 0:
        cemperr += (-0.253 * Eri)**2 * ((0.004 / -0.253)**2 + (Eer / Eri)**2)
    if Bri != 0:
        cemperr += (-0.296 * Bri)**2 * ((0.002 / -0.296)**2 + (Ber / Bri)**2)
    cemperr = cemperr**0.5
    # calculate aggregate UL of empirical system parameters
    ULslvemp = aggregateUL([sempul, aempul, bempul, vempul, lempul, cempul])
    # choose which method to use for system parameters and calculate logKsa, aggregate UL and error
    if type(ULslvqspr) is str or (type(ULslvqspr) is not str and type(ULslvemp) is not str and ULslvqspr < aggregateUL([ULslvemp, ULemptra])):
        if type(ULslvqspr) is str:
            domainnotes.append('system parameters: experimental or user values')
        else:
            domainnotes.append('system parameters: QSPRs applied')
        logKsa = solutedependencies[0]['S'][0] * solventdependencies[0]['s'][0] + \
                 solutedependencies[0]['A'][0] * solventdependencies[0]['a'][0] + \
                 solutedependencies[0]['B'][0] * solventdependencies[0]['b'][0] + \
                 solutedependencies[0]['V'][0] * solventdependencies[0]['v'][0] + \
                 solutedependencies[0]['L'][0] * solventdependencies[0]['l'][0] + \
                 solventdependencies[0]['c'][0]
        logKsaUL = aggregateUL([ULslt, ULslvqspr], roundup=True)
        logKsaerr = (solutedependencies[0]['V'][0]*solventdependencies[0]['v'][2])**2 + \
                    (solutedependencies[0]['L'][0]*solventdependencies[0]['l'][0])**2 * ((solutedependencies[0]['L'][2]/solutedependencies[0]['L'][0])**2 + (solventdependencies[0]['l'][2]/solventdependencies[0]['l'][0])**2) + \
                    solventdependencies[0]['c'][2]**2
        if solutedependencies[0]['S'][0] != 0 and solventdependencies[0]['s'][0] != 0:
            logKsaerr += (solutedependencies[0]['S'][0]*solventdependencies[0]['s'][0])**2 * ((solutedependencies[0]['S'][2]/solutedependencies[0]['S'][0])**2 + (solventdependencies[0]['s'][2]/solventdependencies[0]['s'][0])**2)
        if solutedependencies[0]['A'][0] != 0 and solventdependencies[0]['a'][0] != 0:
            logKsaerr += (solutedependencies[0]['A'][0]*solventdependencies[0]['a'][0])**2 * ((solutedependencies[0]['A'][2]/solutedependencies[0]['A'][0])**2 + (solventdependencies[0]['a'][2]/solventdependencies[0]['a'][0])**2)
        if solutedependencies[0]['B'][0] != 0 and solventdependencies[0]['b'][0] != 0:
            logKsaerr += (solutedependencies[0]['B'][0]*solventdependencies[0]['b'][0])**2 * ((solutedependencies[0]['B'][2]/solutedependencies[0]['B'][0])**2 + (solventdependencies[0]['b'][2]/solventdependencies[0]['b'][0])**2)
        logKsaerr = logKsaerr**0.5
    else:
        if type(ULslvemp) is str:
            domainnotes.append('system parameters: experimental or user solute descriptors plus empirical correlations applied')
            if ULemptra > 1:
                domainnotes.append('solute descriptors out of AD of empirical correlations')
                ULslvemp = aggregateUL([ULslvemp, ULemptra])
        else:
            domainnotes.append('system parameters: solute QSPRs plus empirical correlations applied')
            if ULemptra > 1:
                domainnotes.append('solute descriptors out of AD of empirical correlations')
            ULslvemp = aggregateUL([ULslvemp, ULemptra])
        logKsa = solutedependencies[0]['S'][0] * semp + \
                 solutedependencies[0]['A'][0] * aemp + \
                 solutedependencies[0]['B'][0] * bemp + \
                 solutedependencies[0]['V'][0] * vemp + \
                 solutedependencies[0]['L'][0] * lemp + \
                 cemp
        logKsaUL = aggregateUL([ULslt, ULslvemp], roundup=True)
        logKsaerr = (solutedependencies[0]['V'][0]*vemperr)**2 + \
                    (solutedependencies[0]['L'][0]*lemp)**2 * ((solutedependencies[0]['L'][2]/solutedependencies[0]['L'][0])**2 + (lemperr/lemp)**2) + \
                    cemperr**2
        if solutedependencies[0]['S'][0] != 0 and semp != 0:
            logKsaerr += (solutedependencies[0]['S'][0]*semp)**2 * ((solutedependencies[0]['S'][2]/solutedependencies[0]['S'][0])**2 + (semperr/semp)**2)
        if solutedependencies[0]['A'][0] != 0 and aemp != 0:
            logKsaerr += (solutedependencies[0]['A'][0]*aemp)**2 * ((solutedependencies[0]['A'][2]/solutedependencies[0]['A'][0])**2 + (aemperr/aemp)**2)
        if solutedependencies[0]['B'][0] != 0 and bemp != 0:
            logKsaerr += (solutedependencies[0]['B'][0]*bemp)**2 * ((solutedependencies[0]['B'][2]/solutedependencies[0]['B'][0])**2 + (bemperr/bemp)**2)
        logKsaerr = logKsaerr**0.5
    return [round(logKsa, round_digits)], [logKsaUL], [round(phaseerrorscaling * logKsaerr, round_digits)], '; '.join(domainnotes), citation, units, endpoint

