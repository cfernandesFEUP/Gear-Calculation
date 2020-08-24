## GEAR CALCULATOR ############################################################
import numpy as np,  gearT,  gearM,  maagC,  LStage,  oils,  contact,  heat, \
    bearings, involute
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
## HEAT CONTACT ###############################################################
# alpha_m = 0       #0  -  air
# Re, Pr, hj, hs, hsP, hsW, rc = heat.htc(oil, rohT, niu, k_lub, cp_lub, omega, \
# dsh, ra, rl, alpha_m)
# eqP, eqW, eqPI, eqWI, eqAP, eqAW = heat.heat(rr1, rr2, rB1, rB2, rl, rD1, rD2, ra, avg_qvzp1,\
#                                  avg_qvzp2, qvzp1, qvzp2, a)
## BEARINGS ###################################################################
btype = 'NJ 406'
fab = 0
ngears = 4
Mvl, phi_bl, Msl, Mrr, Mdrag, Grr, Gsl, usl = bearings.pl(btype, frb, fab, n, niu, ubb, ubr)
pvl = ngears*(Mvl[:, :, 0]*omega[0, :] + Mvl[:, :, 1]*omega[1, :])
## TOTAL POWER LOSS (EXCLUDING NO-LOAD) #######################################
pv = pvzp + pvl

## GEOMETRY ###################################################################
size = 200
xP,yP = involute.geo(z[0][0], alpha_t, rb[0][0], ra[0][0], rf[0][0], x[0][0], m, size)
xG,yG = involute.geo(z[0][1], alpha_t, rb[0][1], ra[0][1], rf[0][1], x[0][1], m, size)

import matplotlib.pyplot as plt 
plt.figure(1001)
plt.plot(xP, yP)
# plt.ylim([0,50])
plt.figure(1002)
plt.plot(xG, yG)

## PLOT #######################################################################
# import plot
# plot.fig(xx, vg, qvzp1, qvzp2, avg_qvzp1, avg_qvzp2, lxi, p0, fnx)
# plot.prt(gear, mat, alpha, beta, m, z, x, al, pb, rf, r, rl, ra, epslon_alpha, epslon_beta\
        # , epslon_gama, Req, AB, AC, AD, AE, E, v, cpg, kg, rohg, sigmaHlim, sigmaFlim, Pin, \
        # fbt, fbn, ft, fr, fn, fa, fbear, p0, p0p, Tlub, miu, niu, piezo, HVL, COF, pvzp, pvl,pv)

    
    
    
## CALCULIX ###################################################################


# ds=dsh/1000.0

# import sys
# #sys.path.insert(0, "/home/cfernandes/Gmsh454/lib")
# # MAC
# sys.path.insert(0, "/usr/local/gmsh/lib")
# import gmsh

# model = gmsh.model
# geog = model.geo

# gmsh.initialize()
# gmsh.option.setNumber("General.Terminal", 1)

# model.add("C14")

# # the target mesh size close to the point
# lc = 1e-2

# ## PINION #####################################################################
# h = 10
# R1= 30
# import math
# pointsIPr = []
# pointsIPl = []
# for j in th:
#     pointsIPr.append(geog.addPoint(-rb[0][0]*(np.sin(j)-j*np.cos(j)), rb[0][0]*(np.cos(j)+j*np.sin(j)), 0, h))
# for j in range(len(th)-1, -1, -1):
#     j= th[j]
#     pointsIPl.append(geog.addPoint(rb[0][0]*(np.sin(j-psi_b)-j*np.cos(j-psi_b)), rb[0][0]*(np.cos(j-psi_b)+j*np.sin(j-psi_b)), 0, h))


# linesIPr = geog.addSpline(range(1,len(th)), 1)
# linesIPl = geog.addSpline(range(len(th)+1,2*len(th)), 2)

# centerP = model.geo.addPoint(0,0,0, h)
# # # linesIPl = []
# # # ()
# # # for j in range(len(th)):    
# # #     linesIPl.append(geog.addLine(pointsIPl[j], pointsIPl[(j+1)%len(th)]))
# daP = geog.addCircleArc(pointsIPr[-1], centerP, pointsIPl[0])
# lr = geog.addLine(centerP, pointsIPr[0])
# ll = geog.addLine(pointsIPl[-1], centerP)
# # # raP = model.occ.addCircle(0, 0, 0, ra[0][0])
# # # rloop = geog.addLine(pointsIPr[0], centerP)
# # # lloop = geog.addLine(pointsIPl[0], centerP)
# # # con = geog.addLine(pointsIPr[0], pointsIPl[0])
# curveloop = geog.addCurveLoop([lr,linesIPr, daP, linesIPl, ll], 1)
# # curveloop2 = geog.addCurveLoop([linesIPl], 2)
# # curveloop3 = geog.addCurveLoop([daP], 3)
# # curveloop4 = geog.addCurveLoop([lr], 4)
# # curveloop5 = geog.addCurveLoop([ll], 5)
# disk = geog.addPlaneSurface([curveloop])
# ## WHEEL ######################################################################
# # centerW = model.geo.addPoint(0.0915,0,0, lc, 10)
# # for j in range(2):
# #   points.append(factory.addPoint(ds*math.cos(2*math.pi*j/z[0][0]), ds*math.sin(2*math.pi*j/z[0][0]), 0, h))


# # Create Point for the center of the circle
# # center = model.geo.addPoint(0,0,0, h, 10)
# # Create 3 Points on the circle
# # points = []
# # for j in range(3):
# #   points.append(model.geo.addPoint(R1*math.cos(2*math.pi*j/3), R1*math.sin(2*math.pi*j/3), 0, h))
# # # Create 3 circle arc
# # lines = []
# # for j in range(3):
# #   lines.append(model.geo.addCircleArc(points[j],center,points[(j+1)%3]))
# # # Curveloop and Surface
# # curveloop = model.geo.addCurveLoop([1,2,3])
# # disk = model.geo.addPlaneSurface([curveloop])


# geog.synchronize()
# model.mesh.generate(2)
# gmsh.write("C14.msh")
# gmsh.finalize()