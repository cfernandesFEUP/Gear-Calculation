## GEAR CALCULATOR ############################################################
import numpy as np, gearT, gearM, maagC, LStage, oils, contact, heat
def idt(n,torque):
	## GEAR SELECTION ##
	gear='C14'              # 'C40', '501', '701', '951', 'TPA'
	mat=['POM','POM']       # 'PEEK', 'PA66', 'STEEL' (20MnCr5), 'ADI'
	## TYPE OF GEAR ###############################################################
	alpha,beta,m,z,x,b,dsh,Ra,Rq=gearT.gtype(gear)
	## MAAG CALCULATION ##
	mt,pt,pb,pbt,betab,al,r,rl,ra,rb,rf,alpha_t,alpha_tw,epslon_alpha,epslon_a,\
	epslon_beta,epslon_gama,galpha,galphai,Req,u,T1T2,T1A,T2A,AB,AC,AD,AE,\
	rA1,rA2,rB1,rB2,rD1,rD2=maagC.maag(alpha,beta,m,z,x,b)
	## LINES OF CONTACT ###########################################################
	lxi,xx,rr1,rr2=contact.lines(betab,epslon_alpha,epslon_beta,epslon_gama,rb,T1A,T2A,AE)
	## OPERATING CONDITIONS #######################################################
	Tbulk=50.
	NL=1e6
	n=n/u                 # rpm
	arm='0.35'              # '0.35' or '0.5'
	load=torque#'k02'              # 'k01' up to 'k14 or pinion torque in Nm
	if type(load) is str:
		torque=np.array([[LStage.gtorque(load,arm),u*LStage.gtorque(load,arm)]])
	else:
		torque=np.array([[load,u*load]])
	n=np.array([[u*n,n]])
	omega=np.pi*n/30
	## OIL SELECTION ##############################################################
	oil='dry'
	Tlub=80.0
	Tamb=15.
	if oil=='dry':
		mu=0.28
		rohT,cp_lub,k_lub,beta_lub,piezo,miu,niu,xl=[0.0 for _ in range(8)]
	else:
		mu=0.0
		rohT,cp_lub,k_lub,beta_lub,piezo,miu,niu,xl=oils.astm(oil,Tlub)
	## MATERIAL SELECTION #########################################################
	E,v,cpg,kg,rohg,sigmaHlim,sigmaFlim=gearM.matp(mat,Tbulk,NL)
	## GEAR FORCES ################################################################
	Pin,fbt,fbn,ft,fr,fn,fa,fbear,COF=contact.forces\
	(torque,omega,rb,rl,alpha_tw,betab,Req,Ra,xl,miu,lxi,mu,b)
	## HYERTZ CONTACT #############################################################
	fnx,vt,vri,vr,vg,SRR,Eeff,a,p0,pm,Reff,pvzpx,pvzp,qvzp1,qvzp2,avg_qvzp1,avg_qvzp2,\
	HVL=contact.hertz(lxi,alpha_tw,betab,AE,T1A,T2A,T1T2,rb,E,omega,r,v,fbn,fbt,xx,\
	rr1,Pin,COF,b,pb,kg,cpg,rohg)
	## HEAT CONTACT ###############################################################
	alpha_m=0       #0 - air
	Re, Pr, hj, hs, hsP, hsW, rc = heat.htc(oil, rohT, niu, k_lub, cp_lub, omega, \
    dsh, ra, rl, alpha_m)
    eqP, eqW, eqPI, eqWI, eqAP, eqAW = heat.heat(rr1, rr2, rB1, rB2, rl, rD1, rD2, \
    ra, avg_qvzp1, avg_qvzp2, qvzp1, qvzp2, a)
	## PLOT #######################################################################
	plot.prt(gear,mat,alpha,beta,m,z,x,al,pb,rf,r,rl,ra,epslon_alpha,epslon_beta, \
          epslon_gama,Req,AB,AC,AD,AE,E,v,cpg,kg,rohg,sigmaHlim,sigmaFlim,Pin,\
	        fbt,fbn,ft,fr,fn,fa,fbear,p0,Tlub,miu,niu,piezo,HVL,COF,pvzp)
	return Tamb, rA1, rB1, rl, rD1, eqP, eqW, eqPI, eqWI, eqAP, eqAW