import numpy as np
from scipy import optimize
from GearC import oils
def LCC(ft,z,b,m,x,r,ra,rb,rf,betab,alpha,alpha_t,alpha_tw,u,v,E,rohg,\
        epslon_alpha,epslon_beta,epslon_gama,beta,n,mat,sigmaHlim,\
            sigmaFlim,oil,al,Rz):
    roh40,cp40,k40,beta40,piezo40,miu40,niu40,xl40,ubb40,ubr40=oils.astm(oil,40)
    ## Select the worst operating condition
    ft = ft[-1]
    try:
        n = n[0][-1]
    except:
        n = n[-1]
    KA = 1.00
    zn1 = z[0] / (np.cos(betab) ** 2 * np.cos(beta * np.pi / 180))
    zn2 = z[1] / (np.cos(betab) ** 2 * np.cos(beta * np.pi / 180))
    C1 = 0.04723
    C2 = 0.15551
    C3 = 0.25791
    C4 = -0.00635
    C5 = -0.11654
    C6 = -0.00193
    C7 = -0.24188
    C8 = 0.00529
    C9 = 0.00182
    ql = C1 + C2 / zn1 + C3 / zn2 + C4 * x[0] + (C5 * x[0])\
        / zn1 + C6 * x[1] + (C7 * x[1]) / zn2 + C8 * x[0] ** 2 + C9 * x[1] ** 2
    clth = 1 / ql
    CM = 0.8
    CR = 1
    CB = 1
    if ft * KA / b < 100:
        cl = clth * CM * CR * CB * np.cos(beta * np.pi / 180) * ((ft * KA / b) / 100) **0.25
    else:
        cl = clth * CM * CR * CB * np.cos(beta * np.pi / 180)
    cgama = cl * (0.75 * epslon_alpha + 0.25)
    #print('Tooth stifness: ' + str(cl))
    #print('Mesh stifness: ' + str(cgama))
    ## Internal dynamic factor(KV)
    Ca = 0 ## tip relief

    dm = ra + rf
    #di = dsh
    #q = di/ dm ### rim gears
    ## equivalent mass
    mred = (np.pi / 8) * (dm[0] / (2*rb[0])) ** 2 *\
        (dm[0] ** 2 / (1 / (1e-9 * rohg[0]) + 1 / (1e-9 * rohg[1] * u ** 2))) 
    nE1 = 30000 / (np.pi * z[0]) * np.sqrt(cgama / mred)
    N = n / nE1 ### ressonance ratio
    if ft * KA / b <= 100:
        NS = 0.85
    else:
        NS = 0.5 + 0.35 * np.sqrt(ft * KA / (100 * b))
    if epslon_gama <= 2:
        Cv1 = 0.32 ## pitch deviation effect
        Cv2 = 0.34
        Cv3 = 0.23
        Cv4 = 0.90
        Cv5 = 0.47
        Cv6 = 0.47
    else:
        Cv1 = 0.32 ## pitch deviation effect
        Cv2 = 0.57 / (epslon_gama - 0.3)
        Cv3 = 0.096 / (epslon_gama - 1.56)
        Cv4 = (0.57 - 0.05 * epslon_gama) / (epslon_gama - 1.44)
        Cv5 = 0.47
        Cv6 = 0.12 / (epslon_gama - 1.74)
    if epslon_gama <= 1.5:
        Cv7 = 0.75
    elif epslon_gama > 1.5 or epslon_gama <= 2.5:
        Cv7 = 0.125 * np.sin(np.pi * (epslon_gama - 2)) + 0.875
    elif epslon_gama > 2.5:
        Cv7 = 1
    #Cay1 = (1 / 18) * (sigmaHlim1 / 97 - 18.45) ** 2 + 1.5
    #Cay2 = (1 / 18) * (sigmaHlim2 / 97 - 18.45) ** 2 + 1.5
    #Cay = 0.5 * (Cay1 + Cay2)
    fpb = 0.3 * (m + 0.4 * np.sqrt(2*r[1])) + 4
    falpha = 2.5 * np.sqrt(m) + 0.17 * np.sqrt(2*r[1]) + 0.5
    fbeta = 0.1 * np.sqrt(2*r[1]) + 0.63 * np.sqrt(b) + 4.2
    yp = 0.0758 * fpb
    yf = 0.075 * falpha
    fpb_eff = fpb - yp
    falpha_eff = falpha - yf
    Bp = cl * fpb_eff / (KA * (ft / b))
    Bf = cl * falpha_eff / (KA * (ft / b))
    Bk = abs(1 - cl * Ca / (KA * (ft / b)))
    if N <= NS:
        Kr = (Cv1 * Bp) + (Cv2 * Bf) + (Cv3 * Bk)
        KV = (N * Kr) + 1
    elif N > NS or N <= 1.15:
        KV = (Cv1 * Bp) + (Cv2 * Bf) + (Cv4 * Bk) + 1
    elif N >= 1.5:
        KV = (Cv5 * Bp) + (Cv6 * Bf) + Cv7

    #print('KV:' + str(KV))
    ## Face load factors(KHB and KFB)
    fm = ft * KA + KV
    ## DIN 3990
    ## BHB = 1
    #KlHB = 1.00
    # slHB = 0.05
    # gama = (abs(BHB + KlHB * slHB / (d[0] ** 2) * (d[0] / dsh) ** 4 - 0.3) + 0.3) * (b / d[0]) ** 2
    # fsh0 = 0.023 * gama
    # fsh = fsh0 * fm / b
    fsh = 0
    fma = 0.5 * fbeta
    fbetax = 1.33 * fsh + fma
    ybeta = 0.15 * fbetax
    # xbeta = 0.85
    if ybeta > 6:
        ybeta = 6
    fbetay = fbetax - ybeta
    KHBope = fbetay * cgama / (2 * fm / b)
    if KHBope >= 1:
        KHB = np.sqrt(2 * fbetay * cgama / (fm / b))
        bcalb = np.sqrt((2 * fm / b) / (fbetay * cgama))
    elif KHBope < 1:
        KHB = 1 + fbetay * cgama / (2 * fm / b)
    ## bcalb = 0.5 + (fm / b) / (fbetay * cgama)
    h = np.zeros(2)
    h[0] = ra[0] - rf[0]
    h[1] = ra[1] - rf[1]
    if h[0] > h[1]:
        h = h[0]
    else:
        h = h[1]

    if h / b > (1 / 3):
        hsb = 1 / 3
    else:
        hsb = h / b

    NF = 1 / (1 + (hsb) + (hsb) ** 2)
    KFB = KHB ** NF
    print('KHB: ' + str(KHB))
    print('KFB: ' + str(KFB))
    ## Transverse load factors(KHA and KFA)
    fth = ft * KA * KV * KHB
    yalpha = 0.075 * fpb
    if epslon_gama <= 2:
        KHA = epslon_gama / 2 * (0.9 + 0.4 * cgama * (fpb - yalpha) / (fth / b))
    else:
        KHA = 0.9 + 0.4 * np.sqrt(2 * (epslon_gama - 1) / epslon_gama) * cgama * (fpb - yalpha) / (fth / b)
    if epslon_beta < 1:
        ZEPS = np.sqrt((4 - epslon_alpha) * (1 - epslon_beta) / 3 + epslon_beta / epslon_alpha)
    else:
        ZEPS = np.sqrt(1 / epslon_alpha)
    if KHA > epslon_gama / (epslon_alpha * ZEPS ** 2):
        KHA = epslon_gama / (epslon_alpha * ZEPS ** 2)
    elif KHA <= 1:
        KHA = 1
    KFA = KHA
    print('KHA: ' + str(KHA))
    print('KFA: ' + str(KFA))
    ### Pinion and wheel factors
    M1 = np.tan(alpha_tw) / np.sqrt((np.sqrt(ra[0] ** 2 / rb[0] ** 2 - 1) - 2 * np.pi / z[0])
                              * (np.sqrt(ra[1] ** 2 / rb[1] ** 2 - 1) - (epslon_alpha - 1) * 2 * np.pi / z[1]))
    M2 = np.tan(alpha_tw) / np.sqrt((np.sqrt(ra[1] ** 2 / rb[1] ** 2 - 1) - 2 * np.pi / z[1])
                              * (np.sqrt(ra[0] ** 2 / rb[0] ** 2 - 1) - (epslon_alpha - 1) * 2 * np.pi / z[0]))
    if epslon_beta >= 1:
        ZP = 1
        ZW = 1
    else:
        ZP = M1 - epslon_beta * (M1 - 1)
        ZW = M2 - epslon_beta * (M2 - 1)
        if ZP < 1:
            ZP = 1

        if ZW < 1:
            ZW = 1
    print('ZB: ' + str(ZP))
    print('ZD: ' + str(ZW))
    ## Elasticity factor
    ZE = np.sqrt(1 / (np.pi * (1e6 * (1 - v[0] ** 2) / E[0] + 1e6 * (1 - v[1] ** 2) / E[1])))
    print('ZE: ' + str(ZE))
    ## Zone factor
    ZH = np.sqrt(2 * np.cos(betab) * np.cos(alpha_tw) / (np.cos(alpha_t) ** 2 * np.sin(alpha_tw)))
    print('ZH: ' + str(ZH))
    ## Contact ratio factor
    print('ZEPS: ' + str(ZEPS))
    ## Spiral angle factor
    ZBETA = np.sqrt(np.cos(beta * np.pi / 180))
    print('ZBETA: ' + str(ZBETA))
    ## Lubrication coefficient
    ZL = 0.91 + 0.36 / (1.2 + 134 / niu40) ** 2
    print('ZL: ' + str(ZL))
    vt = np.pi * n * r[0] / 30000
    ZV = 0.93 + 0.14 / (0.8 + 32 / vt) ** 0.5
    #print('ZV: ' + str(ZV))
    ZR = 1.02 * (al ** (1 / 3) / (Rz[0] + Rz[1])) ** 0.08
    print('ZR: ' + str(ZR))
    ## Tooth form factor(YF)
    if beta == 0:
        rfer = 0.375
    else:
        rfer = 0.3
    rfP = rfer * m
    hfP = 1.25 * m
    spr = 0
    zn = np.array([zn1,zn2])
    dn = m*zn
    dbn = dn * np.cos(alpha * np.pi / 180)
    G = rfP / m - hfP / m + x
    Es = np.pi / 4 * m - hfP * np.tan(alpha * np.pi / 180) + \
        spr / np.cos(alpha * np.pi / 180) - (1 - np.sin(alpha * np.pi / 180)) \
            * rfP / np.cos(alpha * np.pi / 180)
    H = 2. / zn * (np.pi / 2 - Es / m) - np.pi / 3
    def eq1(xx):
        return xx - 2 * G[0] / zn[0] * np.tan(xx) + H[0]
    s1 = optimize.brentq(eq1, 0, np.pi / 2.1)
    sol1 = s1
    vu1 = sol1
    def eq2(xx):
        return xx - 2 * G[1] / zn[1] * np.tan(xx) + H[1]
    s2 = optimize.brentq(eq2, 0, np.pi / 2.1)
    sol2 = s2
    vu2 = sol2
    vu = np.array([vu1, vu2])
    sFn = m * (zn * np.sin(np.pi / 3 - vu) + np.sqrt(3) * (G / np.cos(vu) - rfP / m))
    dan = dn + 2*ra - 2*r
    epslon_alphan = epslon_alpha / (np.cos(betab) ** 2)
    FCT1 = ((dan / 2)** 2 - (dbn / 2)** 2) ** 0.5
    FCT2 = np.pi * 2*r * np.cos(beta * np.pi / 180) * np.cos(alpha * np.pi / 180) / z * (epslon_alphan - 1)
    FCT3 = (dbn / 2) ** 2
    den = 2 * ((FCT1 - FCT2) ** 2 + FCT3) ** 0.5
    alpha_en = np.zeros(2)
    gamma_e = np.zeros(2)
    hFe = np.zeros(2)
    Fact1 = np.zeros(2)
    Fact2 = np.zeros(2)
    rF = np.zeros(2)
    qs = np.zeros(2)
    sigmaH = np.zeros(2)
    alpha_en[0] = np.arccos(dbn[0] / den[0])
    alpha_en[1] = np.arccos(dbn[1] / den[1])
    gamma_e[0] = (np.pi / 2 + 2 * x[0] * np.tan(alpha * np.pi / 180)) / zn[0] + np.tan(alpha * np.pi / 180) - alpha * np.pi / 180 - np.tan(
        alpha_en[0]) + alpha_en[0]
    gamma_e[1] = (np.pi / 2 + 2 * x[1] * np.tan(alpha * np.pi / 180)) / zn[1] + np.tan(alpha * np.pi / 180) - alpha * np.pi / 180 - np.tan(
        alpha_en[1]) + alpha_en[1]
    alphaFen = alpha_en - gamma_e
    hFe[0] = 0.5 * m * (
                (np.cos(gamma_e[0]) - np.sin(gamma_e[0]) * np.tan(alphaFen[0])) * den[0] / m - zn[0] * np.cos(np.pi / 3 - vu[0]) - G[0] / np.cos(vu[0]) + rfP / m)
    hFe[1] = 0.5 * m * (
                (np.cos(gamma_e[1]) - np.sin(gamma_e[1]) * np.tan(alphaFen[1])) * den[1] / m - zn[1] * np.cos(np.pi / 3 - vu[1]) - G[1] / np.cos(vu[1]) + rfP / m)
    YF = 6 * (hFe / m) * np.cos(alphaFen) / ((sFn / m) ** 2 * np.cos(alpha * np.pi / 180))
    ## Stress correction factor (YS) ###########################################
    Fact1[0] = 2 * G[0] ** 2
    Fact1[1] = 2 * G[1] ** 2
    Fact2[0] = np.cos(vu[0]) * (zn[0] * np.cos(vu[0]) ** 2 - 2 * G[0])
    Fact2[1] = np.cos(vu[1]) * (zn[1] * np.cos(vu[1]) ** 2 - 2 * G[1])
    rF[0] = m * (rfP / m + Fact1[0] / Fact2[0])
    rF[1] = m * (rfP / m + Fact1[1] / Fact2[1])
    LS = sFn/ hFe
    qs[0] = sFn[0] / (2 * rF[0])
    qs[1] = sFn[1] / (2 * rF[1])
    YS = (1.2 + 0.13 * LS) * qs ** (1. / (1.21 + 2.3 / LS))
    ## Helix angle factor(YB)
    if epslon_beta > 1:
        ebeta = 1
    else:
        ebeta = epslon_beta
    if beta > 30:
        betaDIN = 30
    else:
        betaDIN = beta

    YB = 1 - ebeta * betaDIN / 120
    ## Notch sensitivity factor(YdelT)
    YdelT = 0.9434 + 0.0231 * (1 + 2 * qs) ** 0.5
    # YdelT = 0.44 * YS + 0.12; "static analysis
    #print('YdelT: ' + str(YdelT))
    ## Surface factor(YRrelT)
    YRrelT = 0.957
    # if Rz[0] < 1
    #     RzDIN = 1
    # else
    #     RzDIN = Rz[0]
    # YRrelT = 1.674 - 0.529 * (RzDIN + 1) ** 0.1
    # YRrelT = 5.306 - 4.203 * (RzDIN + 1) ** 0.01
    # YRrelT = 4.299 - 3.259 * (RzDIN + 1) ** 0.0058
    # YRrelT = 1;
    ## static analysis
    #print('YRrelT: ' + str(YRrelT))
    ## PLSTIC AGAINST STEEL ###################################################
    mat = str(mat)
    if mat in 'POM' or mat in 'PA66' or mat in 'PEEK':
        KH = KA
        KF = KA
        ZLUB = ZL * ZV
        SFmin = 2
        YdelT[0] = 1
        YdelT[1] = 1
        YRrelT = 1
    else:
        KH = KA * KV * KHB * KHA
        KF = KA * KV * KFB * KFA
        ZLUB = ZL * ZV * ZR
        SFmin = 1.4
    ## FLANK STRESS ###########################################################
    sigmaH0 = ZE * ZH * ZEPS * ZBETA * np.sqrt(ft * (u + 1) / (b * 2*r[0] * u))
    sigmaH[0] = ZP * sigmaH0 * np.sqrt(KH)
    sigmaH[1] = ZW * sigmaH0 * np.sqrt(KH)
    SHmin = 1.0
    sigmaHP = sigmaHlim / SHmin * ZLUB
    SH = sigmaHP * SHmin / sigmaH
    ## ROOT STRESS ############################################################
    Yst = 2
    sigmaFE = sigmaFlim * Yst
    sigmaF0 = ft / (b * m) * YF * YS * YB
    sigmaF = sigmaF0 * KF
    sigmaFG = sigmaFE * YdelT * YRrelT
    sigmaFP = sigmaFG / SFmin
    SF = np.dot(sigmaF,sigmaFG)/np.dot(sigmaF,sigmaF)
    ## PRINT ##################################################################
    print(sigmaHP)
    print(sigmaH)
    print('Nominal flank pressure: ' + str(sigmaH0))
    print('Flank pressure: ' + str(sigmaH))
    print('Permissible flank pressure: ' + str(sigmaHP))
    print('Flank pressure safety factor (SH): ' + str(SH))
    print('From factor YF: ' + str(YF))
    print('YS: ' + str(YS))
    print('YB: ' + str(YB))
    print('Nominal root stress: ' + str(sigmaF0))
    print('Tooth root stress: ' + str(sigmaF))
    print('Limit tooth root stress: ' + str(sigmaFG))
    print('Permissible root stress: ' + str(sigmaFP))
    print('Root stress safety factor (SF): ' + str(SF))
    return SF, SFmin, SH, SHmin
