## LIBRARY OF MATERIALS #######################################################
def matp(mat, Tbulk, NL):
    import numpy as np
    E, v, cpg, kg, rohg, sigmaHlim, sigmaFlim = [np.zeros(2) for _ in range(7)]
    for i in range(len(E)):
        if mat[i] == 'POM':
            E[i] = 3.2e9        # 2900 MPa (min)  -  3500 MPa (max)
            v[i] = 0.35
            cpg[i] = 1465       # J/kg.K
            kg[i] = 0.3         # W/m.K (0.23 (min) 0.37 (max))
            rohg[i] = 1415      # 1410 (min)  -  1420 (max)
            sigmaHlim[i] = 36 - 0.0012*Tbulk**2 + (1000 - 0.025*Tbulk**2)*NL** - 0.21
            sigmaFlim[i] = 26 - 0.0025*Tbulk**2 + 400*NL** - 0.2
        elif mat[i] == 'PEEK':
            E[i] = 3.65e9
            v[i] = 0.38
            cpg[i] = 1472       # 1443 - 1501
            kg[i] = 0.25        # W/m.K
            rohg[i] = 1320
            sigmaHlim[i] = 36 - 0.0012*Tbulk**2 + (1000 - 0.025*Tbulk**2)*NL** - 0.21 # Nylon (PA66)
            sigmaFlim[i] = 30 - 0.22*Tbulk + (4600 - 900*Tbulk**0.3)*NL**( - 1/3)     # Nylon (PA66)
        elif mat[i] == 'PA66':
            E[i] = 1.85e9       # 1700 MPa (min)  -  2000 MPa (max)
            v[i] = 0.3          # 0.25 - 0.35
            cpg[i] = 1670       # J/kg.K
            kg[i] = 0.26        # W/m.K (0.25 (min) 0.27 (max))
            rohg[i] = 1140      # 1130 (min)  -  1150 (max))
            sigmaHlim[i] = 36 - 0.0012*Tbulk**2 + (1000 - 0.025*Tbulk**2)*NL** - 0.21
            sigmaFlim[i] = 30 - 0.22*Tbulk + (4600 - 900*Tbulk**0.3)*NL**( - 1/3)
        elif mat[i] == 'ADI':
            E[i] = 210e9
            v[i] = 0.26         # 0.22 (min) 0.30 (max)
            cpg[i] = 460.548
            kg[i] = 55          # W/m.K
            rohg[i] = 7850
            sigmaHlim[i] = 1500
            sigmaFlim[i] = 430
        elif mat[i] == 'STEEL':
            E[i] = 206e9
            v[i] = 0.3          # 0.22 (min) 0.30 (max)
            cpg[i] = 465
            kg[i] = 46          # W/m.K
            rohg[i] = 7830
            sigmaHlim[i] = 1500
            sigmaFlim[i] = 430
    return E, v, cpg, kg, rohg, sigmaHlim, sigmaFlim
