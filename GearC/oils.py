## ASTM FORMULAS ##############################################################
import numpy as np
def astm(oil, Tlub):
    base, m_astm, n_astm, alphaT, roh, cp_lub, ubb, ubr = oilt(oil)
    if base == 'MIN':
        s = 0.1390
        t = 0.9904
        xl = 0.846
    elif base == 'PAO':
        s = 0.7382
        t = 0.1335
        xl = 0.666
    elif base == 'PAG':
        s = 0.5489
        t = 0.1485
        xl = 0.585
    elif base == 'EST':
        s = 0.6605
        t = 0.1360
        xl = 0.65
    c_astm = 0.7
    Texp_astm = Tlub + 273.15
    niu = -c_astm+10**(10**(m_astm - n_astm*np.log10(Texp_astm)))
    miu = niu*roh/1000
    piezo = s*(niu**t)*1e-8
    beta_lub = n_astm*(niu + c_astm)*np.log(niu + c_astm)/(niu*Texp_astm)
    k_lub = 120*(1 - 0.0005*(Texp_astm - 273.15)/3)/roh
    rohT = roh+alphaT*roh*(Texp_astm - 15) 
    return rohT, cp_lub, k_lub, beta_lub, piezo, miu, niu, xl, ubb, ubr
## LIBRARY OF LUBRICANTS ######################################################
def oilt(oil):
    if oil == 'MINR':
        base = 'MIN'
        m_astm = 9.0658
        n_astm = 3.4730
        alphaT = -5.8e-4
        roh = 902
        cp_lub = 2306.93
        ubr = (0.0354, 0.0177)
        ubb = (0.0580, 0.056)
    elif oil == 'MINE':
        base = 'PAO'
        m_astm = 7.0478
        n_astm = 2.6635
        alphaT = -6.97e-4
        roh = 893
        cp_lub = 2306.93
        ubr = (0.0435, 0.0078)
        ubb = (0.0440, 0.027)
    elif oil == 'PAOR':
        base = 'PAO'
        m_astm = 7.3514
        n_astm = 2.7865
        alphaT = -5.5e-4
        roh = 859
        cp_lub = 2306.93
        ubr = (0.039, 0.0100)
        ubb = (0.049, 0.044)
    elif oil == 'ESTF':
        base = 'EST'
        m_astm = 7.2610
        n_astm = 2.7493
        alphaT = -6.7e-4
        roh = 957
        cp_lub = 2306.93
    elif oil == 'ESTR':
        base = 'EST'           
        m_astm = 7.5823
        n_astm = 2.8802
        alphaT = -8.1e-4
        roh = 915
        cp_lub = 2306.93
        ubr = (0.0432, 0.0101)
        ubb = (0.060, 0.043)
    elif oil == 'PAGD':
        base = 'PAG'
        m_astm = 5.7597
        n_astm = 2.1512
        alphaT = -7.1e-4
        roh = 1059
        cp_lub = 2306.93
        ubr = (0.0253, 0.0099)
        ubb = (0.0543, 0.044)
    elif oil == 'P150':
        base = 'PAO'
        m_astm = 7.646383
        n_astm = 2.928541
        alphaT = -5.5e-4
        roh = 849
        cp_lub = 2306.93
        ubr = (0.03673, 0.01544)
        ubb = (0.049, 0.044)
    elif oil == 'FVA3':
        base = 'MIN'
        m_astm = 9.2384
        n_astm = 3.5810
        alphaT = -5.5e-4
        roh = 900
        cp_lub = 2000
    elif oil == 'dry':
        base = 'PAO'
        m_astm = 7.646383
        n_astm = 2.928541
        alphaT = -5.5e-4
        roh = 849
        cp_lub = 2306.93
        ubr = (0.03673, 0.01544)
        ubb = (0.049, 0.044)
    return base, m_astm, n_astm, alphaT, roh, cp_lub, ubb, ubr
