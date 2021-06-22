import numpy as np
def field(fbn,Req,Eeff,v,b,rl,alpha_tw,COF):
    xa = np.linspace(-1.5,1.5,200)
    za = np.linspace(0.00001,2.,200)
    
    FNb = 1000*fbn[-1]/b
    RX = Req/1000
    aH = np.sqrt(FNb*RX/(np.pi*Eeff))
    p0 = np.sqrt(FNb*Eeff/(np.pi*RX))
    xH = xa*aH
    zH = za*aH
    
    B = (1/(2*rl[0]/1000*np.sin(alpha_tw))+1/(2*rl[1]/1000*np.sin(alpha_tw)))
    cof = 2*COF.max()
    
    SigmaX, SigmaY, SigmaZ, TauXZ = [np.zeros((len(xH),len(zH))) for _ in range(4)]
    for i in range(len(xH)):
        for j in range(len(zH)):
            
            M = np.sqrt((aH + xH[i])**2 + zH[j]**2)
            N = np.sqrt((aH - xH[i])**2 + zH[j]**2)
            
            phi1 = np.pi*(M + N)/(M*N*np.sqrt(2*M*N + 2*xH[i]**2 + 2*zH[j]**2 -2*aH**2))
            phi2 = np.pi*(M - N)/(M*N*np.sqrt(2*M*N + 2*xH[i]**2 + 2*zH[j]**2 -2*aH**2))
        
            SigmaX[i,j] = -aH*B*Eeff*(zH[j]*((aH**2 + 2*zH[j]**2 + 2*xH[i]**2)*phi1/aH \
                            - 2*np.pi/aH - 3*xH[i]*phi2) + cof*((2*xH[i]**2 - 2*aH**2\
                            - 2*zH[j]**2)*phi2 + 2*np.pi*xH[i]/aH + 2*xH[i]*(aH**2 - xH[i]**2 -zH[j]**2)*phi1/aH))/np.pi
            SigmaY[i,j] = -2*aH*B*Eeff*v[0]*(zH[j]*((aH**2 + zH[j]**2 + xH[i]**2)*phi1/aH \
                            - np.pi/aH - 2*xH[i]*phi2) + cof*((xH[i]**2 - aH**2\
                            - zH[j]**2)*phi2 + np.pi*xH[i]/aH + xH[i]*(aH**2 - xH[i]**2 - zH[j]**2)*phi1/aH))/np.pi
            SigmaZ[i,j] = -aH*B*Eeff*(zH[j]*(aH*phi1 - xH[i]*phi2) + cof*zH[j]**2*phi2)/np.pi
            
            TauXZ[i,j] = -aH*B*Eeff*(zH[j]**2*phi2 + cof*((2*xH[i]**2 + aH**2 + 2*zH[j]**2)*phi1*zH[j]/aH\
                            - 2*np.pi*zH[j]/aH -3*xH[i]*zH[j]*phi2))/np.pi

    Tmax = 0.5*(SigmaX-SigmaZ)
    Toct = np.sqrt((SigmaX-SigmaY)**2 + (SigmaY-SigmaZ)**2 + (SigmaZ-SigmaX)**2)/3
    SvonMises = np.sqrt((SigmaX-SigmaY)**2 + (SigmaY-SigmaZ)**2 + (SigmaZ-SigmaX)**2 + 6*TauXZ**2)/np.sqrt(2)

    import matplotlib.pyplot as plt     
    cmap = plt.get_cmap('jet', 21)
    nc = 21
    
    plt.figure()
    # plt.subplot(2,3,1)
    plt.title(r'$\sigma_{xx}$ / $p_0$')
    cmin = SigmaX.min()/p0
    cmax = SigmaX.max()/p0
    plt.contourf(xH/aH,zH/aH,SigmaX.T/p0,levels=np.linspace(cmin,cmax,nc),cmap=cmap)
    plt.xlabel('x / b')
    plt.ylabel('z / b')
    plt.grid()
    plt.colorbar()
    
    # plt.subplot(2,3,2)
    plt.figure()
    plt.title(r'$\sigma_{yy}$ / $p_0$')
    cmin = SigmaY.min()/p0
    cmax = SigmaY.max()/p0
    plt.contourf(xH/aH,zH/aH,SigmaY.T/p0,levels=np.linspace(cmin,cmax,nc),cmap=cmap)
    plt.xlabel('x / b')
    plt.ylabel('z / b')
    plt.grid()
    plt.colorbar()
    
    plt.figure()
    # plt.subplot(2,3,3)
    plt.title(r'$\sigma_{zz}$ / $p_0$')
    cmin = SigmaZ.min()/p0
    cmax = SigmaZ.max()/p0
    plt.contourf(xH/aH,zH/aH,SigmaZ.T/p0,levels=np.linspace(cmin,cmax,nc),cmap=cmap)
    plt.xlabel('x / b')
    plt.ylabel('z / b')
    plt.grid()
    plt.colorbar()
    
    plt.figure()
    # plt.subplot(2,3,4)
    plt.title(r'$\tau_{xz}$ / $p_0$')
    cmin = TauXZ.min()/p0
    cmax = TauXZ.max()/p0
    plt.contourf(xH/aH,zH/aH,TauXZ.T/p0,levels=np.linspace(cmin,cmax,nc),cmap=cmap)
    plt.xlabel('x / b')
    plt.ylabel('z / b')
    plt.grid()
    plt.colorbar()
    
    plt.figure()
    # plt.subplot(2,3,5)
    plt.title(r'$\tau_{max}$ / $p_0$')
    cmin = Tmax.min()/p0
    cmax = Tmax.max()/p0
    plt.contourf(xH/aH,zH/aH,Tmax.T/p0,levels=np.linspace(cmin,cmax,nc),cmap=cmap)
    plt.xlabel('x / b')
    plt.ylabel('z / b')
    plt.grid()
    plt.colorbar()
    
    plt.figure()
    # plt.subplot(2,3,6)
    plt.title(r'$\tau_{oct}$ / $p_0$')
    cmin = Toct.min()/p0
    cmax = Toct.max()/p0
    plt.contourf(xH/aH,zH/aH,Toct.T/p0,levels=np.linspace(cmin,cmax,nc),cmap=cmap)
    plt.xlabel('x / b')
    plt.ylabel('z / b')
    plt.grid()
    plt.colorbar()
    import matplotlib.ticker as tick
    plt.figure()
    plt.title(r'$\sigma_{von~Mises}$ / $p_0$')
    plt.contourf(xH/aH,zH/aH,SvonMises.T/p0,cmap=cmap,levels = np.linspace(0., 0.32, 15))
    plt.xlabel('x / b')
    plt.ylabel('z / b')
    plt.grid()
    plt.colorbar(format=tick.FormatStrFormatter('%.2f'))
    plt.savefig('logo.png')
    return SigmaX, SigmaY, SigmaZ, TauXZ, Tmax, Toct, SvonMises