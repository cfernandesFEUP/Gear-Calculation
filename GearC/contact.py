import numpy as np
import time
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
def lines(size, b, pbt, betab, epslon_alpha, epslon_beta, epslon_gama, rb, T1A, T2A, AE):
    COND_A = b*np.tan(betab)
    COND_B = epslon_alpha*pbt
    COND_C = epslon_beta*pbt
    COND_D = (epslon_alpha+epslon_beta)*pbt
    ## X DISCRETIZATION
    xx = np.linspace(0., COND_D, size)
    ## FACE WIDTH DISCRETIZATION
    bpos = np.linspace(0, b, len(xx))
    ## NUMBER OF TOOTH FOR THE ANALYSIS
    N = int(np.ceil(epslon_alpha+epslon_beta))
    k = np.arange(-N,N+1)
    ## CREATE VECTORS
    XC = np.zeros([len(k),len(xx),len(bpos)])
    Lh = np.zeros([len(k),len(xx),len(bpos)])
    L = np.zeros([len(xx),len(bpos)])
    ## POSITIONS ALONG PLANE OF ACTION
    t = time.time()
    for i in range(len(k)):
        XC[i] = np.meshgrid(np.add(xx,k[i]*pbt),bpos*np.tan(betab))[0]
    print('TIME LOOP1', time.time() - t)
    ## CONDITIONS FOR THE CALCULATION OF LINES IN CONTACT
    t = time.time()
    if epslon_beta < 1:
        print('TIME LOOP2', time.time() - t)
        indA = np.where((XC >= 0)*(XC < COND_A))
        indB = np.where((XC >= COND_A)*(XC < COND_B))
        indC = np.where((XC >= COND_B)*(XC < COND_D))
        Lh[indA] = XC[indA]/np.sin(betab)
        Lh[indB] = b/np.cos(betab)
        Lh[indC] = b/np.cos(betab) - (XC[indC] - COND_B)/np.sin(betab)
    else:
        indA = np.where((XC >= 0)*(XC < COND_B))
        indB = np.where((XC >= COND_B)*(XC < COND_C))
        indC = np.where((XC >= COND_C)*(XC < COND_D))
        Lh[indA] = XC[indA]/np.sin(betab)
        Lh[indB] = COND_B/np.sin(betab)
        Lh[indC] = COND_B/np.sin(betab) - (XC[indC] - COND_C)/np.sin(betab)
    print('TIME LOOP2', time.time() - t)
    ## CUT THE ARRAYS OVER THE PATH OF CONTACT
    L = np.sum(Lh,axis=0)
    C1 = xx < COND_B
    x_f = xx[C1]/COND_B
    lsum = L[:,C1].T
    lxi = lsum[:,0]/b
    rr1 = np.sqrt((xx*AE + T1A)**2 + rb[0]**2)
    rr2 = np.sqrt((T2A - xx*AE)**2 + rb[1]**2)
    return Lh, L, lxi, lsum, x_f, bpos, rr1, rr2
## HERTZIAN CONTACT ###########################################################
def hertz(lxi, lsum, bpos, alpha_tw, betab, AE, T1A, T2A, T1T2, rb, E, omega, r, v, \
          fbn, fbt, xx, rr1, Pin, COF, b, pbt, kg, cpg, rohg, Req):
    Eeff = 1/((1 - v[0]**2)/E[0] + (1 - v[1]**2)/E[1])
    R1 = (T1A+xx*AE)/1000
    R2 = (T2A-xx*AE)/1000
    vtb = omega[0]*rb[0]/1000
    vt = np.outer(omega[0],np.sqrt((rb[0]/1000)**2 + R1**2))
    vri = np.array([np.outer(omega[0],R1), np.outer(omega[1],R2)])/np.cos(betab)
    R1 = np.tile(R1,(len(bpos),1)).T
    R2 = np.tile(R2,(len(bpos),1)).T
    Reff = 1/((1/R1)+(1/R2))/np.cos(betab)
    fnx, a, p0, pm  = [np.zeros([len(fbn),len(xx),len(bpos)]) for _ in range(4)]
    for i in range(len(fbn)):
        fnx[i] = 1e3*fbn[i]/lsum
        a[i] = np.sqrt((1/np.pi)*fnx[i,:,:]*Reff/Eeff)
        p0[i] = np.sqrt((1/np.pi)*fnx[i]*(Eeff/Reff))
        pm[i] = fnx[i]/(2*a[i])
    p0p = np.array([np.sqrt((2/np.pi)*max(fnx[i,:,0])*(Eeff/(2*Req/1000))) for i in range(len(fbn))])
    term = kg*cpg*rohg
    beta = np.array([term[0]*vri[0],term[1]*vri[1]])
    vr = (vri[0] + vri[1])/2
    vg = abs(vri[0] - vri[1])
    gs1 = abs(vri[0] - vri[1])/vri[0]
    gs2 = abs(vri[0] - vri[1])/vri[1]
    SRR = vg/vr
    fHV = np.zeros([len(bpos),len(xx)])
    for i in range(len(bpos)):
        fHV[i] = fnx[0,:,i]*vg[0]/(fbt[0]*vtb[0])
    HVL = np.trapz(np.trapz(fHV,xx*AE/1000),bpos)/pbt
    pvzp = Pin*HVL*COF.transpose()
    gama = 0.95
    bk1 = beta[0]/(beta[0] + beta[1])
    bk2 = beta[1]/(beta[0] + beta[1])
    pvzpx, fax, qvzp1, qvzp2, avg_qvzp1, avg_qvzp2 = [np.zeros((len(fbn),len(omega[0]),len(xx))) for _ in range(6)]
    for i in range(len(fbn)):
        for j in range(len(omega[0])):
            pvzpx[i,j] = fnx[i,:,0]*vg[j]*COF[i,j]
            fax[i,j] = fnx[i,:,0]*COF[i,j]
            qvzp1[i,j] = gama*bk1[i]*pm[i,:,0]*COF[i,j]*vg[j]
            qvzp2[i,j] = gama*bk2[i]*pm[i,:,0]*COF[i,j]*vg[j]
            avg_qvzp1[i,j] = qvzp1[i,j]*a[i,:,0]*omega[0,j]/(np.pi*vri[0,j])
            avg_qvzp2[i,j] = qvzp2[i,j]*a[i,:,0]*omega[1,j]/(np.pi*vri[1,j])
    return fnx, vt, vri, vr, vg, SRR, Eeff, a, p0, p0p, pm, Reff, R1, R2, pvzpx,\
        fax, pvzp, qvzp1, qvzp2, avg_qvzp1, avg_qvzp2, HVL, bk1, bk2, gs1, gs2