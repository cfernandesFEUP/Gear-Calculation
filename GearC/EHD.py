
def COF(n,vsumc,epslon_alpha,rhoc,fbn,vr,ve,sized,R,young_eq,torque,miu,Blub,piezo,Klub,linesum,TExp,Ra,Rz,Rq):
    Matrix initialization
    h0=cell(sized,1);
    hm=cell(sized,1);
    L=cell(sized,1);
    phiT=cell(sized,1);
    h0T=cell(sized,1);
    hmT=cell(sized,1);
    h0_t=zeros(sized,1);
    hm_t=zeros(sized,1);
    h0_wt=zeros(sized,1);
    hm_wt=zeros(sized,1);
    Lambda_Grubin_wt=zeros(sized,1);
    Lambda_Grubin=zeros(sized,1);
    Lambda_Dowson_wt=zeros(sized,1);
    Sp=zeros(sized,1);
    vrol=zeros(sized,1);
    doleschel=zeros(sized,1);
    xi_doleschel=zeros(sized,1);
    D=zeros(sized,1);
    matsumoto=zeros(sized,1);
    xi_matsumoto=zeros(sized,1);
    fernandes=zeros(sized,1);
    schlenk=zeros(sized,1);
    %% Matsumoto CoF measured
    mu_bl=0.09;
    mu_ehd=0.0093;
    %% Sum velocity on the pitch point
    line=min(min(linesum))/1000;
    R_pp=rhoc/1000;
    for k=1:sized
        vrol(k,1)=mean(vr{k,1}(1,:));
    %% Specific film thickness
        L{k,1}=(Blub(1,k)*(miu(1,k)/1000)*(vr{k,1}).^2)/Klub(1,k);
        phiT{k,1}=(1+0.1*(1+14.8*(ve{k,1}.^.83)).*(L{k,1}.^.64)).^-1;
        h0{k,1}=0.975*R{k,1}.^0.364.*(vrol(k,1).*miu(1,k).*piezo(1,k)).^0.727.*(line*young_eq)^0.091.*fbn(k,1).^-0.091;
        hm{k,1}=1.186*R{k,1}.^0.43.*piezo(1,k).^0.54.*(vrol(k,1).*miu(1,k)).^0.7.*(line)^0.13.*fbn(k,1).^-0.13.*young_eq^-0.03;
        h0T{k,1}=phiT{k,1}.*h0{k,1}*1e6;  
        hmT{k,1}=phiT{k,1}.*hm{k,1}*1e6;
        tmp=abs(R{k,1}(1,:)-2*R_pp);
        [~,pos]=min(tmp);
        h0_t(k,1)=h0T{k,1}(1,pos);
        hm_t(k,1)=hmT{k,1}(1,pos);
        %% Film thickness without thermal correction
        h0_wt(k,1)=mean(mean(h0{k,1}*1e6));
        Lambda_Grubin_wt(k,1)=h0_t(k,1)/(0.5*(Ra(1)+Ra(2)));
        hm_wt(k,1)=mean(mean(hm{k,1}*1e6));
        Lambda_Dowson_wt(k,1)=hm_t(k,1)/(0.5*(Ra(1)+Ra(2)));
        Ra_avg=(Ra(1)+Ra(2))/2;
        %%
        ## GEARS PARAMETERS
        Sp(k,1)=(vsumc(k,1).*(miu(1,k)).*piezo(1,k)^0.5)./fbn(k,1)^0.5;
        Sg=(Ra_avg*1e-6)^2/(epslon_alpha*line*R_pp);
        ## COF FORMULAS
        fernandes(k,1)=0.014*(1/Sp(k,1)).^0.25.*Sg^0.125;
        schlenk(k,1)=0.048*((fbn(k,1)/(line*1000))./(vsumc(k,1).*rhoc)).^0.2.*(1000*miu(1,k))^(-0.05)*Ra_avg^0.25;
        Lambda_Grubin(k,1)=h0_t(k,1)/(0.5*(Ra(1)+Ra(2)));
            if Lambda_Grubin(k,1)>=2:
                xi_doleschel(k,1)=0
            else
                xi_doleschel(k,1)=(1-0.5*Lambda_Grubin(k,1))^
        doleschel(k,1)=xi_doleschel(k,1)*mu_bl+(1-xi_doleschel(k,1))*mu_ehd;
        print('xi from Doleschel', xi_doleschel)
        ## Matsumoto
        D(k,1)=(Rz(1)+Rz(2))/hm_t(k,1);
        xi_matsumoto(k,1)=0.5*log10(D(k,1));
        print('D from Matsumoto', xi_matsumoto)
        matsumoto(k,1)=(1-xi_matsumoto(k,1))*mu_ehd+xi_matsumoto(k,1)*mu_bl;
# speed
# figure(1)
# h12=plot(n(2:end,1),schlenk(2:end,1),'r*');
# set(h12,'linewidth',1)
# hold on
# h13=plot(n(2:end,1),doleschel(2:end,1),'bo');
# set(h13,'linewidth',1)
# hold on
# h15=plot(n(2:end,1),matsumoto(2:end,1),'kx');
# set(h15,'linewidth',2)
# hold on
# % h16=plot(n(2:end,1),fernandes(2:end,1),'r');
# % set(h16,'linewidth',2)
# grid on
# xlabel('n_1 [rpm]')
# ylabel('\mu_{mZ} [-]')
# axis([0 10000 0 0.1])
# set(gcf, 'units', 'centimeters', 'pos', [0 0 8.8 8.8])
# str = 'T_1=133 N\cdotm and \theta_{oil}=40 ^\circC';
# text(1800,0.005,str)
# set(gca,'YTick',[0:0.01:0.1])
# % set(gca,'XTick',[0:400:1200])
# legend('Schlenck (1994)','Doleschel (2002)','Matsumoto (2014)',5);%,'Fernandes'
# %% load
# figure(2)
# h12=plot(torque(2:end,1),schlenk(2:end,1),'r*');
# set(h12,'linewidth',1)
# hold on
# h13=plot(torque(2:end,1),doleschel(2:end,1),'bo');
# set(h13,'linewidth',1)
# hold on
# h15=plot(torque(2:end,1),matsumoto(2:end,1),'kx');
# set(h15,'linewidth',2)
# hold on
# % h16=plot(torque(2:end,1),fernandes(2:end,1),'r');
# % set(h16,'linewidth',2)
# grid on
# xlabel('T_1 [N\cdotm]')
# ylabel('\mu_{mZ} [-]')
# axis([0 300 0 0.1])
# set(gcf, 'units', 'centimeters', 'pos', [0 0 8.8 8.8])
# str = 'n_1=1500 rpm and \theta_{oil}=40 ^\circC';
# text(50,0.005,str)
# set(gca,'YTick',[0:0.01:0.1])
# % set(gca,'XTick',[0:400:1200])
# legend('Schlenck (1994)','Doleschel (2002)','Matsumoto (2014)',5);%,'Fernandes'
# %% temperature
# figure(3)
# h12=plot(TExp(1,:),schlenk,'r*');
# set(h12,'linewidth',1)
# hold on
# h13=plot(TExp(1,:),doleschel,'bo');
# set(h13,'linewidth',1)
# hold on
# h15=plot(TExp(1,:),matsumoto,'kx');
# set(h15,'linewidth',2)
# hold on
# % h16=plot(TExp(1,:),fernandes,'r');
# % set(h16,'linewidth',2)
# grid on
# xlabel('\theta_{oil} [^\circC]')
# ylabel('\mu_{mZ} [-]')
# axis([40 100 0 0.1])
# set(gcf, 'units', 'centimeters', 'pos', [0 0 8.8 8.8])
# str = 'T_1=133 N\cdotm and n_1=1500 rpm';
# text(50,0.005,str)
# set(gca,'YTick',[0:0.01:0.1])
# % set(gca,'XTick',[0:400:1200])
# legend('Schlenck (1994)','Doleschel (2002)','Matsumoto (2014)',5);%,'Fernandes'
# %% film thickness
# % %h0_t=linspace(0,19.7,10000);
# % Lam=h0_t/(0.5*(Ra(1)+Ra(2)));
# % 
# % %hm_t=linspace(0,14.4,10000);
# % 
# % 
# % Dm=(Rz(1)+Rz(2))./(hm_t);
# % 
# % xi_mu=0.5*log10(Dm);
# % 
# % xi_do=0.25*Lam.^2-Lam+1;
# % 
# % xi_do(Lam>=2)=0;

# figure(4)
# h1=semilogx(Lambda_Grubin_wt,xi_doleschel,'-.','linewidth',1.5);
# set(h1,'color','b')
# hold on
# h2=semilogx(Lambda_Dowson_wt,xi_matsumoto,'linewidth',1.5);
# set(h2,'color','k')
# hold on
# xlabel('log \Lambda  [-]')
# grid on
# ylabel('\xi [-]')
# legend('Doleschel (2002)','Matsumoto (2014)',2);
# axis([0.01 2 0 1])
# set(gcf, 'units', 'centimeters', 'pos', [0 0 1.2*8.8 8.8])

# % figure(5)
# % h1=loglog(Lambda_Grubin_wt,Sp,'-.','linewidth',1.5);
# % set(h1,'color','b')
# % hold on
# % h2=loglog(Lambda_Dowson_wt,Sp,'linewidth',1.5);
# % set(h2,'color','k')
# % hold on
# % xlabel('log \Lambda  [-]')
# % grid on
# % ylabel('log Sp [-]')
# % legend('Grubin','Dowson',2);
# % axis([0.01 10 1e-9 1e-6])
# % set(gcf, 'units', 'centimeters', 'pos', [0 0 1.2*8.8 8.8])

# figure(5)
# h1=loglog(Lambda_Grubin_wt,Sp,'-','linewidth',1.5);
# set(h1,'color','b')
# xlabel('log \Lambda  [-]')
# grid on
# hold on
# ylabel('log n [-]')
# legend('F_{bn}=1971','F_{bn}=3922','F_{bn}=5912',2);
# axis([0.01 10 1e-9 1e-6])
# set(gcf, 'units', 'centimeters', 'pos', [0 0 1.2*8.8 8.8])

# figure(6)
# h1=loglog(Lambda_Dowson_wt,Sp,'-','linewidth',1.5);
# set(h1,'color','k')
# xlabel('log \Lambda  [-]')
# grid on
# hold on
# ylabel('log n [-]')
# legend('F_{bn}=1971','F_{bn}=3922','F_{bn}=5912',2);
# axis([0.01 10 1e-9 1e-6])
# set(gcf, 'units', 'centimeters', 'pos', [0 0 1.2*8.8 8.8])
# end

