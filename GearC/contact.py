import numpy as np
## GEAR FORCES AND POWER LOSS #################################################
def forces(torque, omega, rb, rl, alpha_tw, betab, Req, Ra, xl, miu, lxi, mu, b):
    Pin = np.outer(omega[0],torque[0])
    fbt = 1000*torque[0]/rb[0]
    ft = 1000*torque[0]/rl[0]
    fr = fbt*np.sin(alpha_tw)
    fbn = fbt/np.cos(betab)
    fn = ft/np.cos(betab)
    fa = fbt*np.tan(betab);
    fbear = np.sqrt(fr**2 + ft**2)
    frb = fbear/2
    vsumc = 2*omega[0]*rl[0]/1000*np.sin(alpha_tw)
    Ram = (Ra[0] + Ra[1])/2
    lmin = min(lxi)*b
    if mu == 0:
        COF = 0.048*(np.outer(fbn,1/vsumc)/(lmin*Req))**0.2*miu**(-0.05)*Ram**0.25*xl
    else:
        COF = mu*np.ones((fbn.size, omega[0].size))
    return Pin, fbt, fbn, ft, fr, fn, fa, fbear, frb, COF
## LINES OF CONTACT ###########################################################
def lines(betab, epslon_alpha, epslon_beta, epslon_gama, rb, T1A, T2A, AE):
    xx, lxi = [np.linspace(0., 1., 10) for _ in range(2)]
    nn = np.arange(-np.floor(epslon_gama), np.floor(epslon_gama) + 1., 1.)
    xi, xi1, xi2, xi3, li, xxi = [np.zeros((len(nn), len(xx))) for _ in range(6)]
    for i in range(len(nn)):
        xi[i] = xx[:] + nn[i]/epslon_alpha
        xi1[i] = 1*((xi[i] - epslon_gama/epslon_alpha) <= 0)*(xi[i] >= 0)
        xi2[i] = 1*((xi[i] - epslon_beta/epslon_alpha) >= 0)
        xi3[i] = 1*((xi[i] - 1) >= 0)
        if betab>0.0: 
            li[i] = xi1[i]*(xi[i] - xi2[i]*(xi[i] - epslon_beta/epslon_alpha)\
                            - xi3[i]*(xi[i] - 1))/np.sin(betab)
        else:
            li[i] = xi1[i]
    lxi = li.sum(0)
    rr1 = np.sqrt((xx*AE + T1A)**2 + rb[0]**2)
    rr2 = np.sqrt((T2A - xx*AE)**2 + rb[1]**2)
    return lxi, xx, rr1, rr2
## HERTZIAN CONTACT ###########################################################
def hertz(lxi, alpha_tw, betab, AE, T1A, T2A, T1T2, rb, E, omega, r, v, fbn, \
          fbt, xx, rr1, Pin, COF, b, pb, kg, cpg, rohg, Req):
    fnx = np.outer(1e3*fbn, 1/(lxi*b))
    R = np.array([(T1A + xx*AE),(T2A - xx*AE)])*1e-3
    Eeff = 1/((1 - v[0]**2)/E[0] + (1 - v[1]**2)/E[1])
    Reff = 2/((1/R[0]) + (1/R[1]))/np.cos(betab)
    a = np.sqrt((2/np.pi)*fnx*Reff/Eeff)
    p0 = np.sqrt((2/np.pi)*fnx*(Eeff/Reff))
    vt = np.outer(omega[0],np.sqrt((rb[0]/1000)**2 + R[0]**2))
    vri = np.array([np.outer(omega[0],R[0]), np.outer(omega[1],R[1])])/np.cos(betab)
    term = kg*cpg*rohg
    beta = np.array([np.outer(term[0],vri[0]), np.outer(term[1],vri[1])])
    vtb = omega[0]*rb[0]/1000
    p0p = [np.sqrt((2/np.pi)*max(fnx[i])*(Eeff/(2*Req/1000))) for i in range(len(fbn))]
    vr = (vri[0] + vri[1])/2
    vg = abs(vri[1] - vri[0])
    SRR = vg/vr
    pm = fnx/(2*a)
    HVL = np.trapz(fnx[0]*vg[0]/(fbt[0]*vtb[0]), xx*AE/1000)*b/pb
    pvzp = Pin*HVL*COF.transpose()
    gama = 0.95
    bk1 = beta[0]/(beta[0] + beta[1])
    bk2 = beta[1]/(beta[0] + beta[1])
    pvzpx, qvzp1, qvzp2, avg_qvzp1, avg_qvzp2 = [np.zeros((len(fbn), len(omega[0]), len(xx))) for _ in range(5)]
    for i in range(len(fbn)):
        for j in range(len(omega[0])):
            pvzpx[i,j] = fnx[i]*vg[j]*COF[i,j]
            qvzp1[i,j] = gama*bk1[j]*pm[i]*COF[i,j]*vg[j]
            qvzp2[i,j] = gama*bk2[j]*pm[i]*COF[i,j]*vg[j]
            avg_qvzp1[i,j] = qvzp1[i,j]*a[i]*omega[0,j]/(np.pi*vri[0,j])
            avg_qvzp2[i,j] = qvzp2[i,j]*a[i]*omega[1,j]/(np.pi*vri[1,j])
    return fnx, vt, vri, vr, vg, SRR, Eeff, a, p0, p0p, pm, Reff, pvzpx, pvzp, \
    qvzp1, qvzp2, avg_qvzp1, avg_qvzp2, HVL, bk1, bk2
