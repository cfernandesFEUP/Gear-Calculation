import numpy as np
def beart(btype, frb, fab, n, ubb, ubr):
    if btype == 'NJ 406':
        ubear = ubr
        Dr = 90.
        dr = 30.
        dm = (Dr + dr)/2
        R1 = 1*10**-6
        S1 = 0.16
        S2 = 0.0015
        Grr =R1*(dm**2.41)*(frb**0.31)
        Gsl =S1*(dm**0.9)*fab + S2*dm*frb
        Kz = 5.1
        Kl = 0.65
        Krs = 3*10**-8
        Kroll = Kz*Kl*((Dr+dr)/(Dr-dr))*10**-12
        Vm = 0.000149
        B = 23.
        Kdrag = 10*Vm*Kroll*B*(dm**4)
        return dr, Dr, dm, Grr, Gsl, Kz, Kl, Krs, Kdrag, ubear
    elif btype == 'QJ 308N2MA':
        ubear = ubb
        Dr = 90.
        dr = 40.
        dm = (Dr + dr)/2
        R1 = 4.78*10**-7
        R2 = 2.42
        R3 = 1.40*10**-12
        S1 = 1.20*10**-2
        S2 = 0.9
        S3 = 1.40*10**-12
        fgr = R3*dm**4*n**2
        fgs = S3*dm**4*n**2
        Grr = R1*(dm**1.97)*(frb + fgr + R2*fab)**0.54
        Gsl = S1*(dm**0.26)*((frb + fgs)**(4/3) + S2*fab**(4/3))
        Kz = 3.1
        Kl = 0.
        Krs = 3*10**-8
        Kroll = Kz*((Dr + dr)/(Dr - dr))*10**-12
        Vm = 0.000149
        Kdrag = Vm*Kroll*(dm**5)
        return dr, Dr, dm, Grr, Gsl, Kz, Kl, Krs, Kdrag, ubear
def pl(btype, frb, fab, n, niu, ubb, ubr):
    dr, Dr, dm, Grr, Gsl, Kz, Kl, Krs, Kdrag, ubear = beart(btype, frb, fab, n, ubb, ubr)
    phi_bl = 1/(np.exp(1)**((2.6e-8*(n*niu)**1.4)*dm))
    phi_ish = 1/(1 + ((1.84e-9)*((n*dm)**1.28)*(niu**0.64)))
    phi_rs = 1/(np.exp(1)**(Krs*niu*n*(dr + Dr)*(np.sqrt(Kz/(2*(Dr - dr))))))
    usl = (phi_bl*ubear[0])+(1-phi_bl)*ubear[1] # sliding COF
    Msl, Mrr, Mdrag, Mvl = [np.zeros((len(frb), len(n[0, :]), 2)) for _ in range(4)]
    for j in range(len(n[0,:])):
        for i in range(len(frb)):
            Mdrag[i, j, :] = 2*(Kdrag*(n[:, j]**2)) # drag losses
            Msl[i, j, :] = usl[:, j]*Gsl[i] # sliding torque
            Mrr[i, j, :] = phi_ish[:, j]*phi_rs[:, j]*Grr[i]*(niu*n[:, j])**0.6 # rolling torque
    Mvl = (Mrr + Msl + Mdrag)*1e-3 #Total Torque
    return Mvl, phi_bl, Msl, Mrr, Mdrag, Grr, Gsl, usl