import sys
sys.dont_write_bytecode = True
import numpy as np
from GearC import gears,MAAG,contact,LoadStage,oils,material,ISO6336,bearings,plot
## GEAR SELECTION ##################################################################
gear = 'C14'                    # 'C40',  '501',  '701',  '951',  'TPA'
mat = ['STEEL', 'STEEL']        # 'PEEK',  'PA66',  'STEEL' (20MnCr5),  'ADI'
## GEAR FINISHING ##################################################################
Ra = np.array([0.6, 0.6])
Rq = np.array([0.7, 0.7])
Rz = np.array([4.8, 4.8])
## TYPE OF GEAR ####################################################################
alpha, beta, m, z, x, b, dsh = gears.gtype(gear)
## MAAG CALCULATION ##
mt, pt, pb, pbt, betab, al, r, rl, ra, rb, rf, alpha_t, alpha_tw, epslon_alpha,\
epslon_a, epslon_beta, epslon_gama, galpha, galphai, Req, u, T1T2, T1A, T2A, \
AB, AC, AD, AE, rA1, rA2, rB1, rB2, rD1, rD2 = MAAG.calc(alpha, beta, m, z, x, b)
## LINES OF CONTACT ################################################################
size = 100
lxi, lsum, li, xx, bpos, rr1, rr2 = contact.lines(size, b, pbt, betab, \
                epslon_alpha, epslon_beta, epslon_gama, rb, T1A, T2A, AE)
## OPERATING CONDITIONS ############################################################
Tbulk = 50.
NL = 1e6
nmotor = np.array([200., 350., 700., 1050., 1500./u, 1750])# rpm 
load = ['k01','k03','k07','k09'] # 'k01' up to 'k14 or pinion torque in Nm
arm = '0.35'# '0.35' or '0.5' FZG Load Stages
if type(load[0]) is str:
    torqueP = np.array([LoadStage.gtorque(i, arm) for i in load]) 
else:
    torqueP = np.array(load)
torque = np.array([torqueP, u*torqueP])
n = np.array([u*nmotor, nmotor])
omega = np.pi*n/30
## OIL SELECTION ###################################################################
oil = 'P150'
Tlub = 80.0
Tamb = 15.
if oil == 'dry':
    mu = 0.28
    rohT, cp_lub, k_lub, beta_lub, piezo, miu, niu, xl, ubb, ubr = oils.astm(oil, Tlub)
else:
    mu = 0.0
    rohT, cp_lub, k_lub, beta_lub, piezo, miu, niu, xl, ubb, ubr = oils.astm(oil, Tlub)
## MATERIAL SELECTION ##############################################################
E, v, cpg, kg, rohg, sigmaHlim, sigmaFlim = material.matp(mat, Tbulk, NL)
## GEAR FORCES #####################################################################
Pin, fbt, fbn, ft, fr, fn, fa, fbear, frb, COF = contact.forces\
(torque, omega, rb, rl, alpha_tw, betab, Req, Ra, xl, miu, lxi, mu, b)
## HERTZ CONTACT ###################################################################
fnx, vt, vri, vr, vg, SRR, Eeff, a, p0, p0p, pm, Reff, pvzpx, fax, pvzp, qvzp1, qvzp2, \
avg_qvzp1, avg_qvzp2, HVL, bk1, bk2, gs1, gs2 = contact.hertz(lxi, lsum, bpos, alpha_tw, \
betab, AE, T1A, T2A, T1T2, rb, E, omega, r, v, fbn, fbt, xx, rr1, Pin, COF, b, pbt, \
    kg, cpg, rohg, Req)
## BEARINGS ########################################################################
btype = 'NJ 406'
fab = 0
ngears = 4
Mvl, phi_bl, Msl, Mrr, Mdrag, Grr, Gsl, usl = bearings.pl(btype, frb, fab, n, niu, ubb, ubr)
pvl = ngears*(Mvl[:,:,0]*omega[0] + Mvl[:,:,1]*omega[1])
## TOTAL POWER LOSS (EXCLUDING NO-LOAD) ############################################
pv = pvzp + pvl.T
## PLOT ############################################################################
plot.fig(xx, vg, qvzp1, qvzp2, avg_qvzp1, avg_qvzp2, lxi, p0, fnx, load, nmotor)
## ISO 6336/DIN 3990 + VDI 2736 ####################################################
SF, SFmin, SH, SHmin = ISO6336.LCC(ft,z,b,m,x,r,ra,rb,rf,betab,alpha,alpha_t,alpha_tw,u,v,\
E,rohg,epslon_alpha,epslon_beta,epslon_gama,beta,n,mat,sigmaHlim,sigmaFlim,oil,al,Rz)
## PRINT ###########################################################################
print('Gear type:', gear)
print('Gear material:', mat[0], '/', mat[1])
print('Presure angle (\u03B1):', "%.1f" % alpha)
print('Helix angle (\u03B2):',"%.1f" % beta)
print('Module [mm]:',"%.1f" % m)
print('Number of teeth:', int(z[0]),'/', int(z[1]))
print('Profile shift:', x[0],'/', x[1])
print('Axis distance [mm]:',"%.2f" % al)
print('Base pitch [mm]:',"%.2f" % pb)
print('Root radius [mm]:',"%.2f" % rf[0], '/',"%.2f" % rf[1])
print('Reference radius [mm]:',"%.2f" % r[0], '/',"%.2f" % r[1])
print('Pitch radius [mm]:',"%.2f" % rl[0], '/',"%.2f" % rl[1])
print('Tip radius [mm]:', "%.2f" % ra[0], '/',"%.2f" % ra[1])
print('\u03B5_\u03B1:',"%.2f" % epslon_alpha)
print('\u03B5_\u03B2:',"%.2f" % epslon_beta)
print('\u03B5_\u03B3:', "%.2f" % epslon_gama)
print('Equivalent radius [mm]:',"%.2f" % Req)
print('AB [mm]:',"%.2f" % AB)
print('AC [mm]:',"%.2f" % AC)
print('AD [mm]:',"%.2f" % AD)
print('AE [mm]:',"%.2f" % AE)
print('Young Modulus [GPa]:', E[0]/1e9,'/', E[1]/1e9)
print('Poisson Ratio  [-]:', v[0],'/', v[1])
print('Thermal capacity  [J/kg.K]:', cpg[0], '/', cpg[1])
print('Thermal conductivity  [W/m.K]:', kg[0], '/', kg[1])
print('Density  [kg/m3]:', rohg[0], '/', rohg[1])
print('Rooth Strength \u03C3_{Hlim} [MPa]:', sigmaHlim[0],'/', sigmaHlim[1])
print('Flank Strength \u03C3_{Flim} [MPa]:', sigmaFlim[0],'/', sigmaFlim[1])
np.set_printoptions(precision=1)
print('Pin [W]: SPEED x LOAD\n', Pin)
print('Tangential load base circle F_bt [N]:\n', fbt)
print('Tangential load base circle F_bn [N]:\n', fbn)
print('Tangential load F_t [N]:\n', ft)
print('Radial load F_r [N]:\n', fr)
print('Normal load F_n [N]:\n', fn)
print('Axial load F_a [N]:\n', fa)
print('Bearings load (equal distance) F_bearings [N]:\n', fbear)
np.set_printoptions(precision=2)
print('Maximum Hertzian Contact Pressure p_0 [MPa]:\n', \
        np.array([max(p0[0,:,i])/1e6 for i in range(len(fbn))]))
print('Hertzian Contact Pressure Pitch Point p_0 [GPa]:\n', p0p/1e6)
print('Maximum shear stress [MPa]:\n', \
      np.array([max(0.3*p0[0,:,i])/1e6 for i in range(len(fbn))]))
print('Maximum octaedric shear stress [MPa]:\n',\
      np.array([max(0.272*p0[0,:,i])/1e6 for i in range(len(fbn))]))
print('Location of maximum equivalent stress [\u03BCm]:\n', \
      np.array([min(0.7861*a[0,:,i]*1e6) for i in range(len(fbn))]))
print('Oil temperature Tlub [\u00b0C]:',"%.1f" % Tlub)
print('Dynamic Viscosity \u03B7 [mPas]:',"%.2f" % miu)
print('Kinematic Viscosity \u03BD [cSt]:',"%.2f" % niu)
print('Piezoviscosity \u03B1 [1/Pa]:',"%.3f" % piezo)
print('Gear Loss Factor HVL [-]:',"%.4f" % HVL)
np.set_printoptions(precision=4)
print('Gear coefficient of friction \u03BC_mZ [-]:\n', COF)
np.set_printoptions(precision=1)
print('Gear power loss - Pvzp [W]: SPEED x LOAD\n', pvzp)
print('Bearing power loss - Pvl [W]: SPEED x LOAD\n', pvl)
print('Total power loss (excluding no-laod gear losses) Pv [W]: SPEED x LOAD\n', pv)