import numpy as np
def geo(z,alpha_t,rb,ra,rf,x,m,n):
    ## INVOLUTE
    ry = np.linspace(rb, ra, n)
    alphay = np.arccos(rb/ry)
    def inv(alpha):
        return np.tan(alpha) - alpha
    sty = 2*ry*(np.pi/(2*z) + inv(alpha_t) - inv(alphay)) + 2*x*m*np.tan(alpha_t)
    teta1 = sty/(2*ry)
    xI = ry*np.sin(teta1)
    yI = ry*np.cos(teta1)   
    ## ROOT FILLET
    k = np.cos(np.pi/(2*z) + np.tan(alpha_t)-alpha_t - np.tan(alphay) + alphay)/\
        np.sin(np.pi/(2*z) + np.tan(alpha_t)-alpha_t - np.tan(alphay) + alphay)
    teta = np.arctan(-1/k[0])
    rc = (rb-rf)/(1 + np.sin(teta))
    xC = xI[0] + rc*np.cos(teta)
    yC = yI[0] + rc*np.sin(teta)
    alfaC = np.arctan(xC/yC)
    alfaZ = 2*np.pi/z
    tet1 = np.pi/2 + teta
    if alfaZ < alfaC:
        tet0 = alfaC
        tetaF = np.linspace(tet0,tet1,n//2)
        tetaF = np.linspace(tet0,tet1,n//2)
        xF = xC - rc*np.sin(tetaF)
        yF = yC - rc*np.cos(tetaF)
    else:
        tet0 = 0.
        tetaF = np.linspace(tet0,tet1,n//2)
        xF = xC - rc*np.sin(tetaF)
        yF = yC - rc*np.cos(tetaF)
        ## DEDDENDUM CIRCLE
        xd1 = rf*np.sin(np.pi/z)
        xd = np.linspace(xd1, xF[0],n//5)
        yd = np.sqrt(rf**2-xd**2)
        xF = np.concatenate((xd, xF), axis = 0)
        yF = np.concatenate((yd, yF), axis = 0)
    ## ADDENDUM CIRCLE
    xa = np.linspace(xI[-1], -xI[-1],n//10)
    ya = np.sqrt(ra**2-xa**2)
    x = np.concatenate((xF, xI, xa, -xI[::-1], -xF[::-1]), axis=0)
    y = np.concatenate((yF, yI, ya, yI[::-1], yF[::-1]), axis=0)
    return x, y