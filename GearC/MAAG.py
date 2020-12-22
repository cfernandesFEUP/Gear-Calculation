## MAAG GEAR CALCULATION ######################################################
import numpy as np
def calc(alpha, beta, m, z, x, b):
    e_v = 0
    alpha = np.radians(alpha)                               # alpha in rad
    beta = np.radians(beta)                                 # beta in rad
    mt = m/np.cos(beta)                                     # transverse module
    d = z*mt                                                # reference diameter
    alpha_t = np.arctan(np.tan(alpha)/np.cos(beta))         # transverse pressure angle
    betab = np.arcsin(np.sin(beta)*np.cos(alpha))           # base helix angle
    pb = np.pi*m*np.cos(alpha)                              # base pitch
    pt = np.pi*mt                                           # transverse pitch
    pbt = np.pi*m*np.cos(alpha)/np.cos(betab)               # transverse base pitch
    db = d*np.cos(alpha_t)                                  # base diameter
    u = z[1]/z[0]                                           # transmission ratio
    inv_alpha_tw = np.tan(alpha_t) - alpha_t + (2*np.tan(alpha)*np.sum(x)/np.sum(z))
    from scipy import optimize
    def inv(xx):
        return inv_alpha_tw - np.tan(xx) + xx 
    sol  =  optimize.brentq(inv,  0.0,  0.7)
    alpha_tw = sol
    al = np.sum(db)/(2*np.cos(alpha_tw)) + e_v              # working axis distance    
    dl1 = 2*al/(u + 1)                                      # working pitch diameter
    dl2 = 2*u*al/(u + 1)
    dl = np.array([dl1, dl2])
    alpha_tw = np.arccos(np.sum(db)/(2*al))
    haP = 1
    hfP = 1.25
    k = np.sum(z)/2*((((np.tan(alpha_tw) - alpha_tw) - \
                    (np.tan(alpha_t) - alpha_t))/np.tan(alpha)) - \
    1/np.cos(beta)*(np.cos(alpha_t)/np.cos(alpha_tw) - 1))
    da = d + 2*m*(haP + x - k)                              # tip diameter
    df = d + 2*m*(x - hfP)                                  # root diameter
    rb = db/2
    rl = dl/2
    rf = df/2
    ra = da/2
    r = d/2 
    alpha_a = np.arccos(db/da);                             # transverse profile angle at tooth tip
    epslon_a = z*(np.tan(alpha_a) - np.tan(alpha_tw))/(2*np.pi) # addendum contact ratio
    epslon_alpha = np.sum(epslon_a)
    epslon_beta = b*np.tan(betab)/pbt
    epslon_gama = epslon_beta + epslon_alpha
    galphai = rb*(np.tan(alpha_a) - np.tan(alpha_tw))       # length of path of addendum contact individual
    galpha = np.sum(galphai)                                # length of path of contact
    Req = 1/(1/(rl[0]*np.sin(alpha_tw)) + 1/(rl[1]*np.sin(alpha_tw))) # curvature radius on pitch point
    T1T2 = al*np.sin(alpha_tw)
    T1E = np.sqrt(ra[0]**2 - rb[0]**2);
    T2A = np.sqrt(ra[1]**2 - rb[1]**2);
    T1A = T1T2 - T2A
    AE = T1E - T1A
    T1B = T1E - pbt
    T2D = T2A - pbt
    T1D = T1T2 - T2D
    AB = T1B - T1A
    AD = T1D - T1A
    T1C = np.sqrt(rl[0]**2 - rb[0]**2)
    AC = T1C - T1A
    rA1 = np.sqrt(T1A**2 + rb[0]**2)
    rB1 = np.sqrt(T1B**2 + rb[0]**2)
    rD1 = np.sqrt(T1D**2 + rb[0]**2)
    rA2 = np.sqrt((T2A - AE)**2 + rb[1]**2)
    rB2 = np.sqrt((T2A - AD)**2 + rb[1]**2)
    rD2 = np.sqrt((T2A - AB)**2 + rb[1]**2)
    return mt, pt, pb, pbt, betab, al, r, rl, ra, rb, rf, alpha_t, alpha_tw, \
    epslon_alpha, epslon_a, epslon_beta, epslon_gama, galpha, galphai, Req, u,\
    T1T2, T1A, T2A, AB, AC, AD, AE, rA1, rA2, rB1, rB2, rD1, rD2
