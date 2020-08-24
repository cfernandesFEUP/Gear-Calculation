## GEAR CALCULATOR ############################################################
import sys
sys.path.insert(1, '/GearCP/')
import numpy as np,  gearT,  gearM,  maagC,  LStage,  oils,  contact, bearings
## GEAR SELECTION ##
gear = 'C14'                    # 'C40',  '501',  '701',  '951',  'TPA'
mat = ['STEEL', 'STEEL']            # 'PEEK',  'PA66',  'STEEL' (20MnCr5),  'ADI'
## TYPE OF GEAR ###############################################################
alpha, beta, m, z, x, b, dsh, Ra, Rq = gearT.gtype(gear)
## MAAG CALCULATION ##
mt, pt, pb, pbt, betab, al, r, rl, ra, rb, rf, alpha_t, alpha_tw, epslon_alpha,\
epslon_a, epslon_beta, epslon_gama, galpha, galphai, Req, u, T1T2, T1A, T2A, AB,\
 AC, AD, AE, rA1, rA2, rB1, rB2, rD1, rD2 = maagC.maag(alpha, beta, m, z, x, b)
## LINES OF CONTACT ###########################################################
lxi, xx, rr1, rr2 = contact.lines(betab, epslon_alpha, epslon_beta, epslon_gama,\
                                  rb, T1A, T2A, AE)
## OPERATING CONDITIONS #######################################################
Tbulk = 50.
NL = 1e6
nmotor = [200., 350., 700., 1050., 1400., 1850.]# rpm 
arm = '0.35'                    # '0.35' or '0.5'
load = ['k01','k03','k07','k09']# 'k01' up to 'k14 or pinion torque in Nm
torque = np.zeros((len(load), 2))
n = np.zeros((len(nmotor), 2))
for i in range(len(load)):
    if type(load[i]) is str:
        torque[i] = [LStage.gtorque(load[i], arm),\
                                          u*LStage.gtorque(load[i], arm)]
    else:
        torque[i] = [load[i], u*load[i]]
for i in range(len(nmotor)):
    n[i] = [u*nmotor[i], nmotor[i]]
n = n.transpose()
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
E, v, cpg, kg, rohg, sigmaHlim, sigmaFlim = gearM.matp(mat, Tbulk, NL)
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
import plot
plot.fig(xx, vg, qvzp1, qvzp2, avg_qvzp1, avg_qvzp2, lxi, p0, fnx)
plot.prt(gear, mat, alpha, beta, m, z, x, al, pb, rf, r, rl, ra, epslon_alpha, epslon_beta\
       , epslon_gama, Req, AB, AC, AD, AE, E, v, cpg, kg, rohg, sigmaHlim, sigmaFlim, Pin, \
       fbt, fbn, ft, fr, fn, fa, fbear, p0, p0p, Tlub, miu, niu, piezo, HVL, COF, pvzp, pvl,pv)
