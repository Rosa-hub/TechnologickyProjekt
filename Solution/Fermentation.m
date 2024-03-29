function output = Fermentation(m_ext,Yext,mext_out)
clc
close all
VrN=3e+03;
BAF=0;
top=1;
tdisp=0.5;
ybleed=0.1;
S=45;
Vr_Vf=4.5e-5;
tol=1;

while tol>1e-2
    
paramF=fedBatch(S,VrN,Vr_Vf,37,BAF,1);
paramMF=MF(m_ext/paramF.BA);
ro_acid=1/(0.9/1830+0.1/997);
ro_base=880;
n_acid=upravapH(paramMF.Vper,paramF.BA,6,4.5);
n_base=upravapH(mext_out*(1-ybleed),paramF.BA*(1-Yext),4.5,6);
c_acid=0.9*ro_acid/98.08*1000;
c_base=0.28*ro_base/17.031*1000;

V_acid=n_acid/c_acid;
V_base=n_base/c_base;

tR=paramF.tR;

c_BAFMR=paramF.BA*(1-Yext);
c_BA_reg=S/paramF.S*c_BAFMR;
V_FMR=paramF.S*mext_out/S/1000*(1-ybleed);

Vf=m_ext/paramF.BA*(paramF.tR+top);

V_H=Vf/1000-V_FMR;
cBA_Feed=round(V_FMR*c_BA_reg/(V_FMR+V_H),2);

nF=ceil(Vf/(Vf-paramMF.Vper*(paramF.tR-top)));
Vdisp=Vf/nF/tdisp/3600;
VrN=Vf/0.9;
cBA_To_Ext1=round((paramF.BA*paramMF.Vper)/(paramMF.Vper+V_acid)/88.11,2);
Vf_ext=round((m_ext/paramF.BA+V_acid)/1000,2);
tol=abs(cBA_Feed-BAF);
BAF=cBA_Feed;
lifeLK=0.03/paramF.Prod*650*24;

end

output.nR=nF;
output.VR=round(VrN/nF/0.8/1000,1);
output.MFA=paramMF.A;
output.V_ext=Vf_ext;
output.BA_ext=cBA_To_Ext1;
output.Vdisp=Vdisp;
output.Vacid=V_acid;
output.Vbase=V_base;
output.VFeed=V_H/(paramF.tR+top);
output.Prod=paramF.Prod;
output.VLK=paramF.LK/lifeLK;
output.CN=paramF.mN/(paramF.tR+top)/1000;
output.Vrz=Vf/(paramF.tR+top)/1000;
output.BAF=c_BAFMR;
% exportData(paramF,paramMF)


function Param=MF(Vp)
J=(166.4+143.7)/2;

Vper=Vp;
A=Vp/J;

Param.Vper=Vper;
Param.A=round(A,2);



function Param=fedBatch(S,Vr,Vf_Vr,T,BAF,sm)
% 1.1412 C6H12O6	+	0.70964	C5H10O5	+	0.125	NH4NO3	-->	1	CH1.6O0.43N0.25	+	1.5942	C4H8O2	+	0.0368	C2H4O2	+	2.2443	H2	+	2.9455	CO2	+	1.1880	H2O


stech=[-1.1412 -0.7096 -0.1250 1 1.5942 0.0368 2.2443 2.9455 1.1880];
MW=[180.156 150.13 80.043 26.0652 88.11 60.052 2.002 44.01 18.02]; % 8% popol
de=0.701612e-3;
ws=[0.6587 0.3413];
N=20;
kin=[0.48/3600 1.62 372 48.3 5.18];
Ep=0.703;
cLK=75;
VLK=0.1*Vr;
V0=VLK;
cNS=S*ws(1)/MW(1)/stech(1)*stech(3)*MW(3)*1.2;
densFB=1000;


cAi=zeros(1,N);
cAi(N)=S*ws(1);
cBi=zeros(1,N);
cBi(N)=S*ws(2);
cCi=zeros(1,N);
cCi(N)=cNS;
cDi=ones(1,N)*cLK*VLK/(VLK+V0);
cDi(N)=0;
cEi=zeros(1,N);
cEi(N)=BAF;
cFi=zeros(1,N);
cGi=zeros(1,N);
cHi=zeros(1,N);
cIi=zeros(1,N);
Vri=V0;
cIi(N)=0;
Qi=0;

Vf=Vf_Vr*(Vr-Vri);
opt=odeset('Events',@Eventfun);
Feed=[S*ws(1) S*ws(2) cNS 0 BAF 0 0 0 0 0];
IC=[cAi cBi cCi cDi cEi cFi cGi cHi cIi Vri Qi];
tspan=[0 (Vr-Vri)/Vf*200];

[t,c]=ode15s(@kinetModel,tspan,IC,opt,Vf,VLK,cLK,de,N,stech,MW,kin,Feed,IC,Ep,Vr,T);
cA=c(:,1:N);
cB=c(:,N+1:2*N);
cC=c(:,2*N+1:3*N);
cD=c(:,3*N+1:4*N);
cE=c(:,4*N+1:5*N);
cF=c(:,5*N+1:6*N);
cG=c(:,6*N+1:7*N);
cH=c(:,7*N+1:8*N);
cI=c(:,8*N+1:9*N);
Vr=c(:,9*N+1);
Qm=c(:,end);

%Aggitation
dm=0.35;
Np=3;
n=3;
N=250/60;
Qp=Np*N^3*n*dm^5*densFB/1000;

Q=Qm+Qp;
I(1)=0;

for i=2:length(Q)
    I(i)=(Q(i-1)+Q(i))/2*(t(i)-t(i-1));
end

Qs=sum(I)/t(end);

VR=(Vr(end)+VLK)/0.8;

Coilpar=coolingCoil(Qs,T,dm,N,densFB,Vr(end)+VLK);
if sm==0
figure(1)
plot(t/3600,cA(:,end),t/3600,cB(:,end),t/3600,cC(:,end),t/3600,cE(:,end),t/3600,cF(:,end),'LineWidth',2)
legend('Gluc','Xyl','NH4NO3','BA','AA')
title 'Koncentracne profily'
xlabel 't [h]'
ylabel 'Koncentrácia [g/L]'
grid on

figure(2)
plot(t/3600,Q,'LineWidth',2)
title 'Generovanie tepla'
xlabel 't [h]'
ylabel 'Teplo [kW]'
grid on
else 
end

Prod=cE(end,end)/(t(end)/3600+1);
Param.BA=cE(end,end);
Param.S=cA(end,end)+cB(end,end);
Param.VR=VR/1000;
Param.Vf=VR*0.8-VLK;
Param.tR=t(end)/3600;
Param.Vfeed=Vf;
Param.Prod=Prod;
Param.U=Coilpar.U;
Param.L=Coilpar.L;
Param.do=Coilpar.do;
Param.w=Coilpar.w;
Param.nv=Coilpar.nv;
Param.LK=VLK;
Param.mN=cNS*Param.Vf;

% disp(mean(cD(end,:)));

function dxdt=kinetModel(t,x,Vf,VLK,cLK,de,N,stech,MW,kin,Feed,IC,Ep,Vrf,T)
cA=x(1:N);
cB=x(N+1:2*N);
cC=x(2*N+1:3*N);
cD=x(3*N+1:4*N);
cE=x(4*N+1:5*N);
cF=x(5*N+1:6*N);
cG=x(6*N+1:7*N);
cH=x(7*N+1:8*N);
cI=x(8*N+1:9*N);
Vr=x(9*N+1);
Q=x(end);


if Vr>=Vrf
    Vf=0;
end

D=Vf/Vr;

dVr=Vf;

alfa=Vr/VLK;

cS=cA+cB;
cP=cE+cF;
vm=kin(1);
Ks=kin(2);
Ki=kin(3);
Pd=kin(4);
i=kin(5);

for j=1:length(cS)
   
    if cA(j)<1e-8
        cA(j)=1e-8;
    elseif cE(j)<1e-8
        cE(j)=1e-8;
    elseif cB(j)<1e-8
        cB(j)=1e-8;
    elseif cF(j)<1e-8
        cF(j)=1e-8;
    end
    
end

r=vm.*cS./(cS+Ks+(cS.^2./Ki)).*(1-cP./Pd).^i;
R=de/2;
dr=R/(N-1);
Y=stech.*MW/MW(4);
rp=linspace(0,R,N);

C=[cA(N) cB(N) cC(N) cD(N) cE(N) cF(N) cG(N) cH(N)];


DeA=DiffusionCoef(T,1,C,1);
DeB=DiffusionCoef(T,1,C,2);
DeC=DiffusionCoef(T,1,C,3);
DeD=0;
DeE=DiffusionCoef(T,1,C,5);
DeF=DiffusionCoef(T,1,C,6);
DeG=DiffusionCoef(T,1,C,7);
DeH=DiffusionCoef(T,1,C,8);
DeI=0;

De=[DeA DeB DeC DeD DeE DeF DeG DeH DeI];

for j=1:N
    
    if j==1
        
        dAdt(j)=(De(1)*(2*cA(j+1)-2*cA(j))/dr^2+Y(1)*r(j)*cD(j))/Ep;
        dBdt(j)=(De(2)*(2*cB(j+1)-2*cB(j))/dr^2+Y(2)*r(j)*cD(j))/Ep;
        dCdt(j)=(De(3)*(2*cC(j+1)-2*cC(j))/dr^2+Y(3)*r(j)*cD(j))/Ep;
        dDdt(j)=(De(4)*(2*cD(j+1)-2*cD(j))/dr^2+Y(4)*r(j)*cD(j))/Ep;
        dEdt(j)=(De(5)*(2*cE(j+1)-2*cE(j))/dr^2+Y(5)*r(j)*cD(j))/Ep;
        dFdt(j)=(De(6)*(2*cF(j+1)-2*cF(j))/dr^2+Y(6)*r(j)*cD(j))/Ep;
        dGdt(j)=(De(7)*(2*cG(j+1)-2*cG(j))/dr^2+Y(7)*r(j)*cD(j))/Ep;
        dHdt(j)=(De(8)*(2*cH(j+1)-2*cH(j))/dr^2+Y(8)*r(j)*cD(j))/Ep;
        dIdt(j)=(De(9)*(2*cI(j+1)-2*cI(j))/dr^2+Y(9)*r(j)*cD(j))/Ep;
    elseif j==N
        
        dAdt(j)=-3/alfa*De(1)/R*((cA(j)-cA(j-1))/dr)+D*(Feed(1)-cA(j));
        dBdt(j)=-3/alfa*De(2)/R*((cB(j)-cB(j-1))/dr)+D*(Feed(2)-cB(j));
        dCdt(j)=-3/alfa*De(3)/R*((cC(j)-cC(j-1))/dr)+D*(Feed(3)-cC(j));
        dDdt(j)=-3/alfa*De(4)/R*((cD(j)-cD(j-1))/dr)+D*(Feed(4)-cD(j));
        dEdt(j)=-3/alfa*De(5)/R*((cE(j)-cE(j-1))/dr)+D*(Feed(5)-cE(j));
        dFdt(j)=-3/alfa*De(6)/R*((cF(j)-cF(j-1))/dr)+D*(Feed(6)-cF(j));
        dGdt(j)=-3/alfa*De(7)/R*((cG(j)-cG(j-1))/dr)+D*(Feed(7)-cG(j));
        dHdt(j)=-3/alfa*De(8)/R*((cH(j)-cH(j-1))/dr)+D*(Feed(8)-cH(j));
        dIdt(j)=-3/alfa*De(9)/R*((cI(j)-cI(j-1))/dr)+D*(Feed(9)-cI(j));
        

    else
        
        dAdt(j)=(De(1)*((cA(j+1)-2*cA(j)+cA(j-1))/dr^2+2/rp(j)*(cA(j)-cA(j-1))/dr)+Y(1)*r(j)*cD(j))/Ep;
        dBdt(j)=(De(2)*((cB(j+1)-2*cB(j)+cB(j-1))/dr^2+2/rp(j)*(cB(j)-cB(j-1))/dr)+Y(2)*r(j)*cD(j))/Ep;
        dCdt(j)=(De(3)*((cC(j+1)-2*cC(j)+cC(j-1))/dr^2+2/rp(j)*(cC(j)-cC(j-1))/dr)+Y(3)*r(j)*cD(j))/Ep;
        dDdt(j)=(De(4)*((cD(j+1)-2*cD(j)+cD(j-1))/dr^2+2/rp(j)*(cD(j)-cD(j-1))/dr)+Y(4)*r(j)*cD(j))/Ep;
        dEdt(j)=(De(5)*((cE(j+1)-2*cE(j)+cE(j-1))/dr^2+2/rp(j)*(cE(j)-cE(j-1))/dr)+Y(5)*r(j)*cD(j))/Ep;
        dFdt(j)=(De(6)*((cF(j+1)-2*cF(j)+cF(j-1))/dr^2+2/rp(j)*(cF(j)-cF(j-1))/dr)+Y(6)*r(j)*cD(j))/Ep;
        dGdt(j)=(De(7)*((cG(j+1)-2*cG(j)+cG(j-1))/dr^2+2/rp(j)*(cG(j)-cG(j-1))/dr)+Y(7)*r(j)*cD(j))/Ep;
        dHdt(j)=(De(8)*((cH(j+1)-2*cH(j)+cH(j-1))/dr^2+2/rp(j)*(cH(j)-cH(j-1))/dr)+Y(8)*r(j)*cD(j))/Ep;
        dIdt(j)=(De(9)*((cI(j+1)-2*cI(j)+cI(j-1))/dr^2+2/rp(j)*(cI(j)-cI(j-1))/dr)+Y(9)*r(j)*cD(j))/Ep;
        
    end   
    
end


hcA=117.5*(6*4+12-2*6);
hcB=117.5*(5*4+10-5*2);
hcC=117.5*(4-2*2);
hcD=117.5*(4+1.6-0.43*2);
hcE=117.5*(4*4+8-2*2);
hcF=117.5*(2+4-2*2);
hcG=117.5*(2);
hcH=117.5*(0);
hcI=117.5*(0);
hc=[hcA hcB hcC hcD hcE hcF hcG hcH hcI];
dH=sum(stech.*hc);

ksiw=r.*cD;
for i=2:N
  I(i)= (ksiw(i-1)*rp(i-1)^2+ksiw(i)*rp(i)^2)/2*dr; 
end
ksis=3/R^3*sum(I);

dQ=ksis*abs(dH)/alfa;

dxdt=[dAdt dBdt dCdt dDdt dEdt dFdt dGdt dHdt dIdt dVr dQ]';

function param = coolingCoil(Q,T,dm,N,dens,Vr)

tci=15;
tco=T-10;
dtc=tco-tci;
mc=Q/4.2/dtc;
wopt=1.8;
Vc=mc/997;
di=sqrt(4*Vc/pi/wopt)*100/2.54;
BWG=0.259;
do=di+BWG*2;
if mod(do,0.25)<0.25/2
    do_adj=-mod(do,0.25);
    
else
    do_adj=0.25-mod(do,0.25);
end
do=(do+do_adj)*2.54/100;
di=do-2*BWG*2.54/100;
wc=4*Vc/pi/di^2;
Dr=(4*Vr/1000/pi/3)^(1/3);
lam_coil=100;

Rec=wc*di*dens/1e-3;
Prc=6.21;
Nuc=0.023*Rec^0.8*Prc^0.4;
hc=Nuc*0.608/di;

Reh=N*dm^2*1000/10e-3;
Prh=4180*10e-3/0.6;
Nuh=0.187*Reh^0.688*Prh^0.37*(0.69/0.8)^0.11*(Dr/dm)^0.382;
hh=Nuh*0.6/Dr;
LMTD=((T-tci)-(T-tco))/log((T-tci)/(T-tco));
Ul=pi/(1/hc/di+1/2/lam_coil*log(do/di)+1/hh/do);

L=Q*1000/Ul/LMTD;
nv=L/pi/Dr;

Hr=4*Vr/1000/pi/Dr^2;
spc=(Hr-nv*do)/nv;

param.di=di*100/2.54;
param.do=do*100/2.54;
param.w=wc;
param.BWG=BWG;
param.L=round(L,0);
param.nv=ceil(nv);
param.spc=spc;
param.Hr=round(Hr,1);
param.Dr=Dr;
param.U=round(Ul/pi/do,0);

function [val,it,dir]=Eventfun(t,x,Vf,VLK,cLK,de,N,stech,MW,kin,Feed,IC,Ep,Vrf,T)
cA=x(1:N);
cB=x(N+1:2*N);
if x(9*N+1)>=Vrf && (cA(end)+cB(end))<=3
    val=0;
else 
    val=1;
end
it=1;
dir=0;
