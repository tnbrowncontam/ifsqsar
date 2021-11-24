"""Meta QSAR for logKsa"""
import numpy as np
value_names = ('logKsa',)
version = 1
round_digits = 2
units = 'unitless'
components = {'solute': 1, 'solvent': 1}
solute_dependencies_list = ['E', 'S', 'A', 'B', 'V', 'L']
solvent_dependencies_list = ['E', 'S', 'A', 'B', 'V', 'L', 's', 'a', 'b', 'v', 'l', 'c', 'tm']
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

def calculate(solutedependencies, solventdependencies):
    # determine if the solvent is a liquid
    domainnotes = []
    # calculate leverage of the solvent vs. the empirical correlations training dataset
    x = np.array([solventdependencies['E'][0],
                  solventdependencies['S'][0],
                  solventdependencies['A'][0],
                  solventdependencies['B'][0],
                  solventdependencies['V'][0],
                  solventdependencies['L'][0]])
    leverage = np.matmul(np.matmul(x, xtxi), x.T)
    if leverage < 1.5 * 6 / 66:
        ULemptra = 0
    elif leverage < 3 * 6 / 66:
        ULemptra = 1
    elif leverage > 1:
        ULemptra = 3
    else:
        ULemptra = 2
    # calculate MP and MP error
    MP = 53.6 * solventdependencies['E'][0] + \
         27.8 * solventdependencies['S'][0] + \
         107.9 * solventdependencies['A'][0] + \
         -19.2 * solventdependencies['B'][0] + \
         7.1 * solventdependencies['L'][0] + \
         -90.6
    MP = (MP + (solventdependencies['tm'][0] - 273.15)) / 2
    MPerr = (solventdependencies['L'][0] * 7.1)**2 * ((solventdependencies['L'][2] / solventdependencies['L'][0])**2 + (0.7 / 7.1)**2) + \
            3.0**2
    if solventdependencies['E'][0] != 0:
        MPerr += (solventdependencies['E'][0] * 53.6)**2 * ((solventdependencies['E'][2] / solventdependencies['E'][0])**2 + (3.6 / 53.6)**2)
    if solventdependencies['S'][0] != 0:
        MPerr += (solventdependencies['S'][0] * 27.8)**2 * ((solventdependencies['S'][2] / solventdependencies['S'][0])**2 + (4.1 / 27.8)**2)
    if solventdependencies['A'][0] != 0:
        MPerr += (solventdependencies['A'][0] * 107.9)**2 * ((solventdependencies['A'][2] / solventdependencies['A'][0])**2 + (4.9 / 107.9)**2)
    if solventdependencies['B'][0] != 0:
        MPerr += (solventdependencies['B'][0] * -19.2)**2 * ((solventdependencies['B'][2] / solventdependencies['B'][0])**2 + (4.4 / -19.2)**2)
    MPerr = ((0.5 * MPerr**0.5)**2 + (0.5 * solventdependencies['tm'][2])**2)**0.5
    # calculate BP and BPerr
    BP = 13.0 * solventdependencies['E'][0] + \
         43.8 * solventdependencies['S'][0] + \
         59.8 * solventdependencies['A'][0] + \
         18.2 * solventdependencies['B'][0] + \
         26.9 * solventdependencies['V'][0] + \
         29.2 * solventdependencies['L'][0] + \
         -33.0
    BPerr = (solventdependencies['V'][0] * 8.2)**2 + \
            (solventdependencies['L'][0] * 29.2)**2 * ((solventdependencies['L'][2] / solventdependencies['L'][0])**2 + (2.4 / 29.2)**2) + \
            2.8**2
    if solventdependencies['E'][0] != 0:
        BPerr += (solventdependencies['E'][0] * 13.0)**2 * ((solventdependencies['E'][2] / solventdependencies['E'][0])**2 + (3.6 / 13.0)**2)
    if solventdependencies['S'][0] != 0:
        BPerr += (solventdependencies['S'][0] * 43.8)**2 * ((solventdependencies['S'][2] / solventdependencies['S'][0])**2 + (3.3 / 43.8)**2)
    if solventdependencies['A'][0] != 0:
        BPerr += (solventdependencies['A'][0] * 59.8)**2 * ((solventdependencies['A'][2] / solventdependencies['A'][0])**2 + (3.3 / 59.8)**2)
    if solventdependencies['B'][0] != 0:
        BPerr += (solventdependencies['B'][0] * 18.2)**2 * ((solventdependencies['B'][2] / solventdependencies['B'][0])**2 + (2.8 / 18.2)**2)
    BPerr = BPerr**0.5
    # decide if solvent is liquid or not
    if MP+1.96*MPerr <= 25 and BP-1.96*BPerr >= 25 and ULemptra <= 1:
        domainnotes.append('solvent phase: likely liquid, prediction error scaling = 1')
        phaseerrorscaling = 1
    elif (MP <= 25 and BP  >= 25 and ULemptra <= 1) or \
            (MP+1.96*MPerr <= 25 and BP-1.96*BPerr >= 25 and ULemptra > 1) or \
            (MP-1.96*MPerr <= 25 < MP and ULemptra <= 1) or \
            (BP+1.96*BPerr >= 25 > BP and ULemptra <= 1):
        domainnotes.append('solvent phase: maybe liquid, prediction error scaling = 1.25')
        phaseerrorscaling = 1.25
    elif (MP-1.96*MPerr > 25 and ULemptra > 1):
        domainnotes.append('solvent phase: likely solid, prediction error scaling = 2')
        phaseerrorscaling = 2
    elif MP+1.96*MPerr > 25:
        domainnotes.append('solvent phase: maybe solid, prediction error scaling = 1.5')
        phaseerrorscaling = 1.5
    elif (BP+1.96*BPerr < 25 and ULemptra > 1):
        domainnotes.append('solvent phase: likely gas, prediction error scaling = 1.5')
        phaseerrorscaling = 1.5
    elif BP-1.96*BPerr < 25:
        domainnotes.append('solvent phase: maybe gas, prediction error scaling = 1.25')
        phaseerrorscaling = 1.25
    else:
        domainnotes.append('solvent phase: unclassified, prediction error scaling = 2')
        phaseerrorscaling = 2
    # calculate aggregate UL of solute descriptors
    ULslt = 0
    for sltdes in ['S', 'A', 'B', 'L']:
        if solutedependencies[sltdes][1] < 4:
            ULslt += solutedependencies[sltdes][1]**2
        elif solutedependencies[sltdes][1] == 4 and sltdes != 'L':
            ULslt += 1**2
        elif solutedependencies[sltdes][1] == 4 and sltdes == 'L':
            ULslt += 2**2
        elif solutedependencies[sltdes][1] > 4:
            ULslt += 3**2
    ULslt = np.ceil(ULslt/4)**0.5
    # calculate aggregate UL of system parameters from QSPRs
    ULslvqspr = 0
    for slvpar in ['s', 'a', 'b', 'v', 'l', 'c']:
        if solventdependencies[slvpar][1] < 4:
            ULslvqspr += solventdependencies[slvpar][1]**2
        elif solventdependencies[slvpar][1] == 4 and slvpar in ['s', 'a', 'b']:
            ULslvqspr += 1**2
        elif solventdependencies[slvpar][1] == 4 and slvpar == ['v', 'l', 'c']:
            ULslvqspr += 2**2
        elif solventdependencies[slvpar][1] > 4:
            ULslvqspr += 3**2
    ULslvqspr = np.ceil(ULslvqspr/6)**0.5
    # calculate relative solute descriptors, UL and error
    # Vri
    Vri = 1 / solventdependencies['V'][0]
    Vul = 0
    Ver = 0
    # Eri
    Eri = Vri * solventdependencies['E'][0]
    Eul = solventdependencies['E'][1]
    if Eul == 4:
        Eul = 1
    elif Eul in (5, 6):
        Eul = 3
    Eer = Vri * solventdependencies['E'][2]
    # Sri
    Sri = Vri * solventdependencies['S'][0]
    Sul = solventdependencies['S'][1]
    if Sul == 4:
        Sul = 1
    elif Sul in (5, 6):
        Sul = 3
    Ser = Vri * solventdependencies['S'][2]
    # Ari
    Ari = Vri * solventdependencies['A'][0]
    Aul = solventdependencies['A'][1]
    if Aul == 4:
        Aul = 1
    elif Aul in (5, 6):
        Aul = 3
    Aer = Vri * solventdependencies['A'][2]
    # Bri
    Bri = Vri * solventdependencies['B'][0]
    Bul = solventdependencies['B'][1]
    if Bul == 4:
        Bul = 1
    elif Bul in (5, 6):
        Bul = 3
    Ber = Vri * solventdependencies['B'][2]
    # Lri
    Lri = Vri * solventdependencies['L'][0]
    Lul = solventdependencies['L'][1]
    if Lul == 4:
        Lul = 1
    elif Lul in (5, 6):
        Lul = 3
    Ler = Vri * solventdependencies['L'][2]
    # IAB
    IAB = (Ari * Bri)**0.5
    IABul = ((Aul**2 + Bul**2)/2)**0.5
    IABer = 0
    if IAB > 0:
        IABer = Ari * Bri * ((Aer/Ari)**2 + (Ber/Bri)**2)**0.5
        IABer = IAB * 0.5 * IABer / (Ari * Bri)
    # IVL
    IVL = -np.log((np.exp(-Vri) + np.exp(-Lri)) / 2)
    IVLul = ((Vul**2 + Lul**2)/2)**0.5
    IVLer = 0.5 * abs(np.exp(-Lri)) * Ler / (0.5 * (np.exp(-Vri) + np.exp(-Lri)))
    # calculate empirical system parameters, UL and error
    # s emp
    if Eri <= 0 and Sri <= 0:
        semp = 0
        sempul = 1
        semperr = 0
    elif Eri <= 0.15 and Sri <= 0.15:
        semp = 0.16
        sempul = 1
        semperr = semp / 1.96
    else:
        semp = 1.327 * Sri - 0.405 * Vri - 0.239 * Lri - 0.507 * IAB + 1.871
        sempul = ((Sul**2 + Lul**2 + IABul**2) / 3)**0.5
        if semp < 0:
            semp = 0
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
        aempul = 1
        aemperr = 0
    elif Bri < 0.1:
        aemp = 0.54
        aempul = 1
        aemperr = aemp / 1.96
    else:
        aemp = 2.845 * Bri + 0.421 * IAB - 1.582 * Sri + 0.607 * Lri
        aempul = ((Bul**2 + IABul**2 + Sul**2 + Lul**2) / 4)**0.5
        if aemp < 0:
            aemp = 0
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
        bempul = 1
        bemperr = 0
    else:
        bemp = 0.747 * Ari + 0.378 * IAB - 0.228 * Sri + 0.312 * Vri
        bempul = ((Aul**2 + IABul**2 + Sul**2) / 3)**0.5
        bemperr = 0.312**2 * Vri**2 * 0.007**2
        if bemp < 0:
            bemp = 0
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
    vempul = ((Eul**2 + Sul**2 + IVLul**2 + IABul**2) / 4)**0.5
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
    lempul = ((Eul**2 + Sul**2 + IVLul**2 + IABul**2) / 4)**0.5
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
    cempul = ((Eul**2 + Bul**2 + IVLul**2) / 3)**0.5
    cemperr = (0.157 * IVL)**2 * ((0.004 / 0.157)**2 + (IVLer / IVL)**2) + 0.005**2
    if Eri != 0:
        cemperr += (-0.253 * Eri)**2 * ((0.004 / -0.253)**2 + (Eer / Eri)**2)
    if Bri != 0:
        cemperr += (-0.296 * Bri)**2 * ((0.002 / -0.296)**2 + (Ber / Bri)**2)
    cemperr = cemperr**0.5
    # calculate aggregate UL of empirical system parameters
    ULslvemp = ((sempul**2 + aempul**2 + bempul**2 + vempul**2 + lempul**2 + cempul**2)/6)**0.5
    ULslvemp = np.ceil(((ULslvemp**2 + ULemptra**2)/2)**0.5)
    # choose which method to use for system parameters and calculate logKsa, aggregate UL and error
    if ULslvqspr < ULslvemp:
        domainnotes.append('system parameters: QSPRs applied')
        logKsa = solutedependencies['S'][0] * solventdependencies['s'][0] + \
                 solutedependencies['A'][0] * solventdependencies['a'][0] + \
                 solutedependencies['B'][0] * solventdependencies['b'][0] + \
                 solutedependencies['V'][0] * solventdependencies['v'][0] + \
                 solutedependencies['L'][0] * solventdependencies['l'][0] + \
                 solventdependencies['c'][0]
        logKsaUL = int(np.ceil(((ULslt**2 + ULslvqspr**2)/2)**0.5))
        logKsaerr = (solutedependencies['V'][0]*solventdependencies['v'][2])**2 + \
                    (solutedependencies['L'][0]*solventdependencies['l'][0])**2 * ((solutedependencies['L'][2]/solutedependencies['L'][0])**2 + (solventdependencies['l'][2]/solventdependencies['l'][0])**2) + \
                    solventdependencies['c'][2]**2
        if solutedependencies['S'][0] != 0 and solventdependencies['s'][0] != 0:
            logKsaerr += (solutedependencies['S'][0]*solventdependencies['s'][0])**2 * ((solutedependencies['S'][2]/solutedependencies['S'][0])**2 + (solventdependencies['s'][2]/solventdependencies['s'][0])**2)
        if solutedependencies['A'][0] != 0 and solventdependencies['a'][0] != 0:
            logKsaerr += (solutedependencies['A'][0]*solventdependencies['a'][0])**2 * ((solutedependencies['A'][2]/solutedependencies['A'][0])**2 + (solventdependencies['a'][2]/solventdependencies['a'][0])**2)
        if solutedependencies['B'][0] != 0 and solventdependencies['b'][0] != 0:
            logKsaerr += (solutedependencies['B'][0]*solventdependencies['b'][0])**2 * ((solutedependencies['B'][2]/solutedependencies['B'][0])**2 + (solventdependencies['b'][2]/solventdependencies['b'][0])**2)
        logKsaerr = logKsaerr**0.5
    else:
        domainnotes.append('system parameters: solute QSPRs plus empirical correlations applied')
        logKsa = solutedependencies['S'][0] * semp + \
                 solutedependencies['A'][0] * aemp + \
                 solutedependencies['B'][0] * bemp + \
                 solutedependencies['V'][0] * vemp + \
                 solutedependencies['L'][0] * lemp + \
                 cemp
        logKsaUL = int(np.ceil(((ULslt**2 + ULslvemp**2) / 2) ** 0.5))
        logKsaerr = (solutedependencies['V'][0]*vemperr)**2 + \
                    (solutedependencies['L'][0]*lemp)**2 * ((solutedependencies['L'][2]/solutedependencies['L'][0])**2 + (lemperr/lemp)**2) + \
                    cemperr**2
        if solutedependencies['S'][0] != 0 and semp != 0:
            logKsaerr += (solutedependencies['S'][0]*semp)**2 * ((solutedependencies['S'][2]/solutedependencies['S'][0])**2 + (semperr/semp)**2)
        if solutedependencies['A'][0] != 0 and aemp != 0:
            logKsaerr += (solutedependencies['A'][0]*aemp)**2 * ((solutedependencies['A'][2]/solutedependencies['A'][0])**2 + (aemperr/aemp)**2)
        if solutedependencies['B'][0] != 0 and bemp != 0:
            logKsaerr += (solutedependencies['B'][0]*bemp)**2 * ((solutedependencies['B'][2]/solutedependencies['B'][0])**2 + (bemperr/bemp)**2)
        logKsaerr = logKsaerr**0.5
    return round(logKsa, round_digits), logKsaUL, round(phaseerrorscaling * logKsaerr, round_digits), ', '.join(domainnotes)

