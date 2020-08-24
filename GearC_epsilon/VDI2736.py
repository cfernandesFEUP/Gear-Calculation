function [sigmaF,cl]=VDI(ft,b,m,d,u,v,E,rohg,epslon_alpha,...
    epslon_beta,epslon_gama,beta,betab,alpha,alphat,alpha_tw,n,...
    z1,z2,x1,x2,da,db,df,mat,sigmaHlim,sigmaFlim,niu40,a_l,Rz)
%% Application factor (KA) - method B
KA=1.00; %uniform
%% Stifness -method B
zn1=z1/(cos(betab)^2*cos(beta*pi/180));
zn2=z2/(cos(betab)^2*cos(beta*pi/180));
C1=0.04723;
C2=0.15551;
C3=0.25791;
C4=-0.00635;
C5=-0.11654;
C6=-0.00193;
C7=-0.24188;
C8=0.00529;
C9=0.00182;
ql=C1+C2/zn1+C3/zn2+C4*x1+(C5*x1)/zn1+C6*x2+(C7*x2)/zn2+C8*x1^2+C9*x2^2;
clth=1/ql;
CM=0.8;
CR=1;
CB=1;
if ft*KA/b<100
    cl=clth*CM*CR*CB*cos(beta*pi/180)*((ft*KA/b)/100)^0.25;
else
    cl=clth*CM*CR*CB*cos(beta*pi/180);
end
cgama=cl*(0.75*epslon_alpha+0.25);
disp('Tooth stifness')
disp(cl)
disp('Mesh stifness')
disp(cgama)
%% Internal dynamic factor (KV)
Ca=0; %tip relief
dm=0.5*(da+df);
% di=dsh; 
% q=di./dm; %rim gears
mred=(pi/8)*(dm(1)/db(1))^2*(dm(1)^2/(1/(1e-9*rohg(1))+1/(1e-9*rohg(2)*u^2)));  %equivalent mass
nE1=30000/(pi*z1)*sqrt(cgama/mred)
N=n(1)/nE1;     %ressonance ratio
if ft*KA/b<=100
    NS=0.85;
else
    NS=0.5+0.35*sqrt(ft*KA/(100*b));
end
if epslon_gama<=2
    Cv1=0.32;   %pitch deviation effect
    Cv2=0.34;
    Cv3=0.23;
    Cv4=0.90;
    Cv5=0.47;
    Cv6=0.47;
else
    Cv1=0.32;   %pitch deviation effect
    Cv2=0.57/(epslon_gama-0.3);
    Cv3=0.096/(epslon_gama-1.56);
    Cv4=(0.57-0.05*epslon_gama)/(epslon_gama-1.44);
    Cv5=0.47;
    Cv6=0.12/(epslon_gama-1.74);
end
if epslon_gama<=1.5
    Cv7=0.75;  
elseif epslon_gama>1.5 || epslon_gama<=2.5
    Cv7=0.125*sin(pi*(epslon_gama-2))+0.875;
elseif epslon_gama>2.5
    Cv7=1;
end
%Cay1=(1/18)*(sigmaHlim1/97-18.45)^2+1.5;
%Cay2=(1/18)*(sigmaHlim2/97-18.45)^2+1.5;
%Cay=0.5*(Cay1+Cay2);
fpb=0.3*(m+0.4*sqrt(d(2)))+4;
falpha=2.5*sqrt(m)+0.17*sqrt(d(2))+0.5;
fbeta=0.1*sqrt(d(2))+0.63*sqrt(b)+4.2;
yp=0.0758*fpb;
yf=0.075*falpha;
fpb_eff=fpb-yp;
falpha_eff=falpha-yf;
Bp=cl*fpb_eff/(KA*(ft/b));
Bf=cl*falpha_eff/(KA*(ft/b));
Bk=abs(1-cl*Ca/(KA*(ft/b)));
if N<=NS
    Kr=(Cv1*Bp)+(Cv2*Bf)+(Cv3*Bk);
    KV=(N*Kr)+1;
elseif N>NS || N<=1.15
    KV=(Cv1*Bp)+(Cv2*Bf)+(Cv4*Bk)+1;
elseif N>=1.5
    KV=(Cv5*Bp)+(Cv6*Bf)+Cv7;
end
disp('KV');
disp(KV);
%% Face load factors (KHB and KFB)
fm=ft*KA+KV;
% DIN 3990
%BHB=1;
%KlHB=1.00;
%slHB=0.05;
%gama=(abs(BHB+KlHB*slHB/(d(1)^2)*(d(1)/dsh)^4-0.3)+0.3)*(b/d(1))^2;
%fsh0=0.023*gama;
%fsh=fsh0*fm/b;
fsh=0;
fma=0.5*fbeta;
fbetax=1.33*fsh+fma;
ybeta=0.15*fbetax;
%xbeta=0.85;
if ybeta>6
    ybeta=6;
end
fbetay=fbetax-ybeta;
KHBope=fbetay*cgama/(2*fm/b);
if KHBope>=1
    KHB=sqrt(2*fbetay*cgama/(fm/b));
    %bcalb=sqrt((2*fm/b)/(fbetay*cgama));
elseif KHBope<1
    KHB=1+fbetay*cgama/(2*fm/b);
    %bcalb=0.5+(fm/b)/(fbetay*cgama);
end
h(1)=(da(1)-df(1))/2;
h(2)=(da(2)-df(2))/2;
if h(1)>h(2)
    h=h(1);
else
    h=h(2);
end
if h/b>(1/3)
    hsb=1/3;
else
    hsb=h/b;
end
NF=1/(1+(hsb)+(hsb)^2);
KFB=KHB^NF;
disp('KHB');
disp(KHB);
disp('KFB');
disp(KFB);
%% Transverse load factors (KHA and KFA)
fth=ft*KA*KV*KHB;
yalpha=0.075*fpb;
if epslon_gama<=2
    KHA=epslon_gama/2*(0.9+0.4*cgama*(fpb-yalpha)/(fth/b));
else
    KHA=0.9+0.4*sqrt(2*(epslon_gama-1)/epslon_gama)*cgama*(fpb-yalpha)/(fth/b);
end
if epslon_beta<1                                   
    ZEPS=sqrt((4-epslon_alpha)*(1-epslon_beta)/3+epslon_beta/epslon_alpha);
elseif epslon_beta>=1
    ZEPS=sqrt(1/epslon_alpha);
end
if KHA>epslon_gama/(epslon_alpha*ZEPS^2)
    KHA=epslon_gama/(epslon_alpha*ZEPS^2);
elseif KHA<=1
    KHA=1;
end
KFA=KHA;
disp('KHA');
disp(KHA);
disp('KFA');
disp(KFA);
%% Pinion and wheel factors
M1=tan(alpha_tw)/sqrt((sqrt(da(1)^2/db(1)^2-1)-2*pi/z1)...
    *(sqrt(da(2)^2/db(2)^2-1)-(epslon_alpha-1)*2*pi/z2));
M2=tan(alpha_tw)/sqrt((sqrt(da(2)^2/db(2)^2-1)-2*pi/z2)...
    *(sqrt(da(1)^2/db(1)^2-1)-(epslon_alpha-1)*2*pi/z1));
if epslon_beta>=1
    ZP=1;
    ZW=1;
else
    ZP=M1-epslon_beta*(M1-1);
    ZW=M2-epslon_beta*(M2-1);
    if ZP<1
        ZP=1;
    end
    if ZW<1
        ZW=1;
    end
end
disp('ZB');
disp(ZP);
disp('ZD');
disp(ZW);
%% Elasticity factor
ZE=sqrt(1/(pi*(1e6*(1-v(1)^2)/E(1)+1e6*(1-v(2)^2)/E(2))));
disp('ZE');
disp(ZE);
%% Zone factor
ZH=sqrt(2*cos(betab)*cos(alpha_tw)/(cos(alphat)^2*sin(alpha_tw)));
disp('ZH');
disp(ZH);
%% Contact ratio factor
disp('ZEPS');
disp(ZEPS);
%% Spiral angle factor
ZBETA=sqrt(cos(beta*pi/180));                      
disp('ZBETA');
disp(ZBETA);
%% Lubrication coefficient
ZL=0.91+0.36/(1.2+134/niu40)^2;
disp('ZL');
disp(ZL);
vt=pi*n(1)*d(1)/60000;
ZV=0.93+0.14/(0.8+32/vt)^0.5;
disp('ZV');
disp(ZV);
ZR=1.02*(a_l^(1/3)/(Rz(1)+Rz(2)))^0.08;
disp('ZR');
disp(ZR);
%% Tooth form factor (YF)
if beta==0
    rfer=0.375;
else
    rfer=0.3;
end
rfP=rfer*m;
hfP=1.25*m;
spr=0;
z=[z1 z2];
zn=[zn1 zn2];
xs=[x1 x2];
dn=m*zn;
dbn=dn*cos(alpha*pi/180);
G=rfP/m-hfP/m+xs;
Es=pi/4*m-hfP*tan(alpha*pi/180)+spr/cos(alpha*pi/180)-(1-sin(alpha*pi/180))*rfP/cos(alpha*pi/180);
H=2./zn.*(pi/2-Es/m)-pi/3;
syms xx1
eq1=xx1==2*G(1)/zn(1)*tan(xx1)-H(1);
s1 = feval(symengine, 'numeric::realroots', eq1, 'xx1=0..pi/2.1', 1e-8);
sol1=double(s1(1));
vu1=double(sol1(1));
syms xx2
eq2=xx2==2*G(2)/zn(2)*tan(xx2)-H(2);
s2 = feval(symengine, 'numeric::realroots', eq2, 'xx2=0..pi/2.1', 1e-8);
sol2=double(s2(1));
vu2=double(sol2(1));
vu=[vu1 vu2];
sFn=m*(zn.*sin(pi/3-vu)+sqrt(3)*(G./cos(vu)-rfP/m));
dan=dn+da-d;
epslon_alphan=epslon_alpha/(cos(betab)^2);
FCT1=((dan/2).^2-(dbn/2).^2).^0.5;
FCT2=pi*d.*cos(beta*pi/180)*cos(alpha*pi/180)./z.*(epslon_alphan-1);
FCT3=(dbn/2).^2;
den=2*((FCT1-FCT2).^2+FCT3).^0.5;
alpha_en(1)=acos(dbn(1)/den(1));
alpha_en(2)=acos(dbn(2)/den(2));
gamma_e(1)=(pi/2+2*xs(1)*tan(alpha*pi/180))./zn(1)+tan(alpha*pi/180)-alpha*pi/180-tan(alpha_en(1))+alpha_en(1);
gamma_e(2)=(pi/2+2*xs(2)*tan(alpha*pi/180))./zn(2)+tan(alpha*pi/180)-alpha*pi/180-tan(alpha_en(2))+alpha_en(2);
alphaFen=alpha_en-gamma_e;
hFe(1)=0.5*m*((cos(gamma_e(1))-sin(gamma_e(1))*tan(alphaFen(1)))*den(1)/m-zn(1)*cos(pi/3-vu(1))-G(1)/cos(vu(1))+rfP/m);
hFe(2)=0.5*m*((cos(gamma_e(2))-sin(gamma_e(2))*tan(alphaFen(2)))*den(2)/m-zn(2)*cos(pi/3-vu(2))-G(2)/cos(vu(2))+rfP/m);
YF=6*(hFe/m).*cos(alphaFen)./((sFn/m).^2*cos(alpha*pi/180));
disp('YF');
disp(YF);
%% Stress correction factor (YS)
Fact1(1)=2*G(1)^2;
Fact1(2)=2*G(2)^2;
Fact2(1)=cos(vu(1))*(zn(1)*cos(vu(1))^2-2*G(1));
Fact2(2)=cos(vu(2))*(zn(2)*cos(vu(2))^2-2*G(2));
rF(1)=m*(rfP/m+Fact1(1)/Fact2(1));
rF(2)=m*(rfP/m+Fact1(2)/Fact2(2));
LS=sFn./hFe;
qs(1)=sFn(1)/(2*rF(1));
qs(2)=sFn(2)/(2*rF(2));
YS=(1.2+0.13*LS).*qs.^(1./(1.21+2.3./LS));
disp('YS');
disp(YS);
%% Helix angle factor (YB)
if epslon_beta>1
    ebeta=1;
else
    ebeta=epslon_beta;
end
if beta>30
    betaDIN=30;
else
    betaDIN=beta;
end
YB=1-ebeta*betaDIN/120;
disp('YB');
disp(YB);
%% Notch sensitivity factor (YdelT)
YdelT=0.9434+0.0231*(1+2*qs).^0.5;
%YdelT=0.44*YS+0.12; %static analysis
disp('YdelT');
disp(YdelT);
%% Surface factor (YRrelT)
YRrelT=0.957;
% if Rz(1)<1
%     RzDIN=1;
% else
%     RzDIN=Rz(1);
% end
%YRrelT=1.674-0.529*(RzDIN+1)^0.1;
%YRrelT=5.306-4.203*(RzDIN+1)^0.01;
%YRrelT=4.299-3.259*(RzDIN+1)^0.0058;
%YRrelT=1; static analysis
disp('YRrelT');
disp(YRrelT);
%% Plastic and Steel
if strcmp(mat(1),'POM') || strcmp(mat(1),'PA66') || strcmp(mat(1),'PEEK')
    KH=KA;
    KF=KA;
    ZLUB=ZL*ZV;
    SFmin=2;
    YdelT(1)=1;
    YdelT(2)=1;
    YRrelT=1;
else
    KH=KA*KV*KHB*KHA;
    KF=KA*KV*KFB*KFA;
    ZLUB=ZL*ZV*ZR;
    SFmin=1.4;
end
%% Flank pressure
sigmaH0=ZE*ZH*ZEPS*ZBETA*sqrt(ft*(u+1)/(b*d(1)*u));
sigmaH(1)=ZP*sigmaH0*sqrt(KH);
sigmaH(2)=ZW*sigmaH0*sqrt(KH);
SHmin=1.0;
sigmaHP=sigmaHlim/SHmin*ZLUB;
SH=sigmaHP*SHmin./sigmaH;
disp('Nominal flank pressure:');
disp(sigmaH0);
disp('Flank pressure:');
disp(sigmaH);
disp('Permissible flank pressure:');
disp(sigmaHP);
disp('Flank pressure safety factor (SH):');
disp(SH);
%% Root stress
Yst=2;
sigmaFE=sigmaFlim*Yst;
disp('Nominal root stress:');
sigmaF0=ft/(b*m)*YF.*YS.*YB;
disp(sigmaF0);
disp('Tooth root stress:');
sigmaF=sigmaF0.*KF;
disp(sigmaF);
disp('Limit tooth root stress:');
sigmaFG=sigmaFE*YdelT*YRrelT;
disp(sigmaFG);
disp('Permissible root stress:');
sigmaFP=sigmaFG/SFmin;
disp(sigmaFP);
disp('Root stress safety factor (SF):');
SF=sigmaFG/sigmaF;
disp(SF);
end

