## PLOT #######################################################################
import matplotlib.pyplot as plt 
def fig(xx, vg, qvzp1, qvzp2, avg_qvzp1, avg_qvzp2, lxi, p0, fnx, load, nmotor):
    for i in range(len(load)):
        plt.figure(1)
        plt.plot(xx,  1e-3*fnx[0,:,i], label=load[i])
        plt.ylabel(r'$Fn\left(x\right)$ / Nmm$^{-1}$')
        plt.xlabel('\u03B6')
        plt.title('Load per face width')
        plt.grid(True)
        plt.legend()
        plt.figure(2)
        plt.plot(xx, p0[0,:,i]*1e-9,label=load[i])
        plt.ylabel(r'$p_0\left(x\right)$ / GPa')
        plt.title('Hertzian Contact Pressure')
        plt.xlabel('\u03B6')
        plt.grid(True)
        plt.legend()
    for j in range(len(nmotor)):
        plt.figure(3)
        plt.plot(xx, abs(vg[j, :]),label=str(int(nmotor[j]))+' rpm')
        plt.ylabel(r'$v_g\left(x\right)$ / ms$^{-1}$')
        plt.title('Sliding Speed')
        plt.xlabel('\u03B6')
        plt.grid(True)
        plt.legend()
    plt.figure(4)
    plt.plot(xx, lxi, 'k-')
    plt.xlabel('\u03B6')
    plt.ylabel(r'$L_t\left(x\right)/b$')
    plt.grid(True)
    plt.title('Total contact length over face width')
    plt.figure(5)
    plt.plot(xx, fnx[0,:,-1]/max(fnx[0,:,-1]), 'k-')
    plt.xlabel('\u03B6')
    plt.ylabel(r'$\overline{F_N}~\left(x\right)$')
    plt.title('Dimensionless normal load')
    plt.grid(True)
    plt.figure(6)
    plt.plot(xx, p0[0,:,-1]/max(p0[0,:,-1]), 'k-')
    plt.xlabel('\u03B6')
    plt.ylabel(r'$\overline{p_0}~\left(x\right)$')
    plt.title('Dimensionless contact pressure')
    plt.grid(True)
    plt.figure(7)
    MAX = max(max(qvzp1[:,-1,-1]),max(qvzp2[:,-1,-1]))
    plt.plot(xx, qvzp1[:,-1,-1]/MAX, 'k-', label='pinion')
    plt.plot(xx, qvzp2[:,-1,-1]/MAX, 'r--', label='wheel')
    plt.xlabel('\u03B6')
    plt.ylabel(r'$q_{vzp}~\left(x\right)$')
    plt.title('Dimensionless heat Flux')
    plt.legend()
    plt.grid(True)
    plt.figure(8) 
    MAX_AVG = max(max(avg_qvzp1[:,-1,-1]),max(avg_qvzp2[:,-1,-1]))    
    plt.plot(xx, avg_qvzp1[:,-1,-1]/MAX_AVG, 'k-', label='pinion')
    plt.plot(xx, avg_qvzp2[:,-1,-1]/MAX_AVG, 'r--',  label='wheel')
    plt.legend(('Pinion', 'Wheel'))
    plt.xlabel('\u03B6')
    plt.ylabel(r'$\overline{q_{vzp}}~\left(x\right)$')
    plt.grid(True)
    plt.title('Dimensionless average heat flux')
    plt.show()