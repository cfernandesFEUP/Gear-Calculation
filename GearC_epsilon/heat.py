import numpy as np
## HEAT TRANSFER COEFFICIENTS #################################################
def htc(oil, rohT, niu, k_lub, cp_lub, omega, dsh, ra, rl, alpha_m):
    k_air = 2.85e-2
    cp_air = 1009
    roh_air = 1.067
    niu_air = 18.9e-6
    alpha_air = k_air/(roh_air*cp_air)
    Nu, Lc, roh_mist, niu_mist, k_mist, cp_mist, alpha_mist, Re, Pr, d, rc, hj,\
    hs = [np.zeros((1000, 2)) for _ in range(13)]
    for i in range(0, 2):
        rc[:, i] = np.linspace(dsh/2, ra[0, i], 1000)
#        if oil == 'dry':
        Pr[:, i] = niu_air/alpha_air
        Re[:, i] = 4*np.pi*(rc[:, i]/1000)**2*omega[i, 0]/niu_air
        hj[:, i] = 0.228*Re[:, i]**0.731*Pr[:, i]**(1/3)*k_air/(1e-3*rl[0, i])
        Lc[:,i] = k_air*(omega[i, 0]/niu_air)**0.5
        Nu[:, i] = 0.0188*Re[:, i]**0.8 
        hs[:, i] = Nu[:, i]*Lc[:, i] #Hartnett
        hs[:, i] = 0.335*Lc[:, i]
#        hs[:, i] = 0.36*Re[:, i]**0.5*k_mist[:, i]/(2*rc[:, i])
#        else:
#            d[:, i] = rc[:, i]/(1e-3*ra[0, i])
#            roh_mist[:, i] = (1-alpha_m*d[:, i])*roh_air + alpha_m*d[:, i]*rohT
#            niu_mist[:, i] = (1-alpha_m*d[:, i])*niu_air + alpha_m*d[:, i]*niu/1e6
#            k_mist[:, i] = (1-alpha_m*d[:, i])*k_air + alpha_m*d[:, i]*k_lub
#            cp_mist[:, i] = (1-alpha_m*d[:, i])*cp_air + alpha_m*d[:, i]*cp_lub
#            alpha_mist[:, i] = k_mist[:, i]/(roh_mist[:, i]*cp_mist[:, i])
#            Re[:, i] = (2*rc[:, i])**2*omega[0, i]/niu_mist[:, i]
#            Pr[:, i] = niu_mist[:, i]/alpha_mist[:, i]
#            hj[:, i] = 0.228*Re[:, i]**0.731*Pr[:, i]**(1/3)*k_mist[:, i]/(1e-3*rl[0, i])
            
    hsP = np.polyfit(rc[:, 0], hs[:, 0], 4)
    hsW = np.polyfit(rc[:, 1], hs[:, 1], 4)
#    import matplotlib.pyplot as plt
#    plt.plot(rc[:, 0], hs[:, 0])
#    plt.plot(rc[:, 0], np.polyval(hsP, rc[:, 0]))
#    plt.plot(rc[:, 1], hs[:, 1])
#    plt.plot(rc[:, 1], np.polyval(hsW, rc[:, 1]))
    return Re, Pr, hj, hs, hsP, hsW, rc
## HEAT CONTACT ###############################################################
def heat(rr1, rr2, rB1, rB2, rl, rD1, rD2, ra, avg_qvzp1, avg_qvzp2, qvzp1, qvzp2, a):
    xp1 = (rr1[:] < rB1)
    xp2 = (rr1[:] >= rB1)*(rr1[:] < rl[0, 0])
    xp3 = (rr1[:] >= rl[0, 0])*(rr1[:] < rD1)
    xp4 = (rr1[:] >= rD1)*(rr1[:] <= ra[0, 0])
    xw1 = (rr2[:] < rB2)
    xw2 = (rr2[:] >= rB2)*(rr2[:] < rl[0, 1])
    xw3 = (rr2[:] >= rl[0, 1])*(rr2[:] < rD2)
    xw4 = (rr2[:] >= rD2)*(rr2[:] <= ra[0, 1])
    xp=[xp1, xp2, xp3, xp4]
    xw=[xw1, xw2, xw3, xw4]
    eqP, eqW, eqPI, eqWI, eqAP, eqAW = [[] for _ in range(6)]
    for i in range(4):
        ordpp = (i==3)*2+(i!=3)*3
        ordpw = (i==0)*2+(i!=0)*3
        eqP.append(np.polyfit(np.trim_zeros(rr1*xp[i]), np.trim_zeros(xp[i]*avg_qvzp1[0, 0, :]), ordpp))
        eqW.append(np.polyfit(np.trim_zeros(rr2*xw[i]), np.trim_zeros(xw[i]*avg_qvzp2[0, 0, :]), ordpw))
        eqPI.append(np.polyfit(np.trim_zeros(rr1*xp[i]), np.trim_zeros(xp[i]*qvzp1[0, 0, :]), ordpp))
        eqWI.append(np.polyfit(np.trim_zeros(rr2*xw[i]), np.trim_zeros(xw[i]*qvzp2[0, 0, :]), ordpw))
        eqAP.append(np.polyfit(np.trim_zeros(rr1*xp[i]),  np.trim_zeros(xp[i]*a[0, :]), 3))
        eqAW.append(np.polyfit(np.trim_zeros(rr2*xw[i]),  np.trim_zeros(xw[i]*a[0, :]), 3))
    return eqP, eqW, eqPI, eqWI, eqAP, eqAW