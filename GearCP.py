import numpy as np
from GearC import gears, MAAG, contact, LoadStage, oils, material, bearings, plot
## GEAR SELECTION ##
gear = 'C40'                    # 'C40',  '501',  '701',  '951',  'TPA'
mat = ['STEEL', 'STEEL']        # 'PEEK',  'PA66',  'STEEL' (20MnCr5),  'ADI'
## TYPE OF GEAR ###############################################################
alpha, beta, m, z, x, b, dsh, Ra, Rq = gears.gtype(gear)
## MAAG CALCULATION ##
mt, pt, pb, pbt, betab, al, r, rl, ra, rb, rf, alpha_t, alpha_tw, epslon_alpha,\
epslon_a, epslon_beta, epslon_gama, galpha, galphai, Req, u, T1T2, T1A, T2A, \
AB, AC, AD, AE, rA1, rA2, rB1, rB2, rD1, rD2 = MAAG.calc(alpha, beta, m, z, x, b)
## LINES OF CONTACT ###########################################################
lxi, xx, rr1, rr2 = contact.lines(betab, epslon_alpha, epslon_beta, epslon_gama,\
                                    rb, T1A, T2A, AE)
## OPERATING CONDITIONS #######################################################
Tbulk = 50.
NL = 1e6
nmotor = np.array([200., 350., 700., 1050., 1400., 1850.])# rpm 
arm = '0.35'# '0.35' or '0.5' FZG Load Stages
load = ['k01','k03','k07','k09']# 'k01' up to 'k14 or pinion torque in Nm
if type(load[0]) is str:
    torqueP = np.array([LoadStage.gtorque(i, arm) for i in load])
else:
    torqueP = np.array(load)
torque = np.array([torqueP, u*torqueP])
n = np.array([u*nmotor, nmotor])
omega = np.pi*n/30
## OIL SELECTION ##############################################################
oil = 'dry'
Tlub = 80.0
Tamb = 15.
if oil == 'dry':
    mu = 0.28
    rohT, cp_lub, k_lub, beta_lub, piezo, miu, niu, xl, ubb, ubr = oils.astm(oil, Tlub)
else:
    mu = 0.0
    rohT, cp_lub, k_lub, beta_lub, piezo, miu, niu, xl, ubb, ubr = oils.astm(oil, Tlub)
## MATERIAL SELECTION #########################################################
E, v, cpg, kg, rohg, sigmaHlim, sigmaFlim = material.matp(mat, Tbulk, NL)
## GEAR FORCES ################################################################
Pin, fbt, fbn, ft, fr, fn, fa, fbear, frb, COF = contact.forces\
(torque, omega, rb, rl, alpha_tw, betab, Req, Ra, xl, miu, lxi, mu, b)
## HERTZ CONTACT #############################################################
fnx, vt, vri, vr, vg, SRR, Eeff, a, p0, p0p, pm, Reff, pvzpx, pvzp, qvzp1, qvzp2, \
avg_qvzp1, avg_qvzp2, HVL, bk1, bk2 = contact.hertz(lxi, alpha_tw, betab, AE, T1A, T2A, \
T1T2, rb, E, omega, r, v, fbn, fbt, xx, rr1, Pin, COF, b, pb, kg, cpg, rohg, Req)
## BEARINGS ###################################################################
btype = 'NJ 406'
fab = 0
ngears = 4
Mvl, phi_bl, Msl, Mrr, Mdrag, Grr, Gsl, usl = bearings.pl(btype, frb, fab, n, niu, ubb, ubr)
pvl = ngears*(Mvl[:, :, 0]*omega[0, :] + Mvl[:, :, 1]*omega[1, :])
## TOTAL POWER LOSS (EXCLUDING NO-LOAD) #######################################
pv = pvzp + pvl
## PLOT #######################################################################
plot.fig(xx, vg, qvzp1, qvzp2, avg_qvzp1, avg_qvzp2, lxi, p0, fnx)
np.set_printoptions(precision=3)
print('Gear type: ', gear)
print('Gear material: ', mat)
print('Presure angle (\u03B1): ', "%.1f" % alpha)
print('Helix angle (\u03B2): ',"%.1f" % beta)
print('Module [mm]: ',"%.1f" % m)
print('Number of teeth: ', int(z[0]),'/', int(z[1]))
print('Profile shift: ', x[rl0],'/', x[1])
print('Axis distance [mm]: ',"%.2f" % al)
print('Base pitch [mm]: ',"%.2f" % pb)
print('Root radius [mm]: ', rf[0,0],'/', rf[0,1])
print('Reference radius [mm]: ', r[0,0],'/', r[0,1])
print('Pitch radius [mm]: ', rl[0,0],'/', rl[0,1])
print('Tip radius [mm]: ', ra[0,0],'/', ra[0,1])
print('\u03B5_\u03B1: ',"%.2f" % epslon_alpha)
print('\u03B5_\u03B2: ',"%.2f" % epslon_beta)
print('\u03B5_\u03B3: ', "%.2f" %epslon_gama)
print('Equivalent radius [mm]: ',"%.2f" % Req)
print('AB [mm]: ',"%.2f" % AB)
print('AC [mm]: ',"%.2f" % AC)
print('AD [mm]: ',"%.2f" % AD)
print('AE [mm]: ',"%.2f" % AE)
print('Young Modulus [GPa]: ', E[0]/1e9,'/', E[1]/1e9)
print('Poisson Ratio  [-]: ', v[0],'/', v[1])
print('Thermal capacity  [J/kg.K]: ', cpg)
print('Thermal conductivity  [W/m.K]: ', kg)
print('Density  [kg/m3]: ', rohg)
print('\u03C3_{Hlim} [MPa]: ', sigmaHlim)
print('\u03C3_{Flim} [MPa]: ', sigmaFlim)
print('Pin [W]: ', Pin)
print('Tangential load base circle F_bt [N]: ', fbt)
print('Tangential load base circle F_bn [N]: ', fbn)
print('Tangential load F_t [N]: ', ft)
print('Radial load F_r [N]: ', fr)
print('Normal load F_n [N]: ', fn)
print('Axial load F_a [N]: ', fa)
print('Bearings load (equal distance) F_bearings [N]: ', fbear)
# print('p_0 [GPa]: ', max(p0[0])/1e9)
# print('p_0p [GPa]: ', p0p)
print('Tlub [\u03B1 C] ', Tlub)
print('\u03B7 [mPas] ', miu)
print('\u03BD [cSt] ', niu)
print('\u03B1 [1/Pa] ', piezo)
print('HVL [-] ', HVL)
print('\u03BC_mZ [-] ', COF)
print('Pvzp [W] ', pvzp)
print('Pvl [W] ', pvl)
print('Pv [W] ', pv)
