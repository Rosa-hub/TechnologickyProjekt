function Fermentation
clc
close all

param=fedBatch(45,5000,1e-5,37)

function Param=fedBatch(S,Vr,Vf_Vr,T)
% 0.9782	C6H12O6	+	0.6086	C5H10O5	+	0.1250	NH4NO3	-->	1 CH1.6O0.43N0.25	+	1.2450	C4H8O2	+	0.0287	C2H4O2	+	2.7902	H2	+	2.8732	CO2	+	0.5617	H2O

stech=[-0.9872 -0.6086 -0.1250 1 1.2450 0.0287 2.7902 2.8732 0.5617];
MW=[180.156 150.13 80.043 26.0652 88.11 60.052 2.002 44.01 18.02]; % 8% popol
de=0.701612e-3;
ws=[0.6587 0.3413];
N=20;
kin=[0.48/3600 1.62 372 48.3 5.18];
Ep=0.703;
cLK=50;
VLK=0.2*Vr;
V0=0.5*VLK;
Vf=Vf_Vr*(Vr-VLK-V0);
cNS=S*ws(1)/MW(1)/stech(1)*stech(3)*MW(3)*1.2;


cAi=zeros(1,N);
cAi(N)=S*ws(1);
cBi=zeros(1,N);
cBi(N)=S*ws(2);
cCi=zeros(1,N);
cCi(N)=cNS;
cDi=ones(1,N)*cLK;
cDi(N)=0;
cEi=zeros(1,N);
cFi=zeros(1,N);
cGi=zeros(1,N);
cHi=zeros(1,N);
cIi=zeros(1,N);
Vri=VLK+V0;
Qi=0;


De=[1e-12 1e-12 1e-12 0 1e-12 1e-12 1e-12 1e-12 1e-12 0];
Feed=[S*ws(1) S*ws(2) cNS 0 0 0 0 0 0 0];
IC=[cAi cBi cCi cDi cEi cFi cGi cHi cIi Vri Qi];
tspan=[0 (Vr-VLK-V0)/Vf];

[t,c]=ode15s(@kinetModel,tspan,IC,[],Vf,De,VLK,cLK,de,N,stech,MW,kin,Feed,IC,Ep);
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
Q=c(:,end);



plot(t/3600,cA(:,end),t/3600,cB(:,end),t/3600,cC(:,end),t/3600,cE(:,end),t/3600,cF(:,end),'LineWidth',2)
legend('Gluc','Xyl','NH4NO3','BA','AA')

% plot(t/3600,Q)

Param=cE(end,end);


function dxdt=kinetModel(t,x,Vf,De,VLK,cLK,de,N,stech,MW,kin,Feed,IC,Ep)
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
   
    if cS(j)<1e-8
        cS(j)=1e-8;
    elseif cP(j)<1e-8
        cP(j)=1e-8;
    elseif cD(j)<1e-8
        cD(j)=1e-8;
    end
    
end

r=vm.*cS./(cS+Ks+(cS.^2./Ki)).*(1-cP./Pd).^i;
R=de/2;
dr=R/(N-1);
Y=stech.*MW/MW(4);
rp=linspace(0,R,N);

for j=1:N
    
    if j==1
        
        dAdt(j)=De(1)/Ep*((2*cA(j+1)-2*cA(j))/dr^2+Y(1)*r(j)*cD(j));
        dBdt(j)=De(2)/Ep*((2*cB(j+1)-2*cB(j))/dr^2+Y(2)*r(j)*cD(j));
        dCdt(j)=De(3)/Ep*((2*cC(j+1)-2*cC(j))/dr^2+Y(3)*r(j)*cD(j));
        dDdt(j)=De(4)/Ep*((2*cD(j+1)-2*cD(j))/dr^2+r(j)*cD(j));
        dEdt(j)=De(5)/Ep*((2*cE(j+1)-2*cE(j))/dr^2+Y(5)*r(j)*cD(j));
        dFdt(j)=De(6)/Ep*((2*cF(j+1)-2*cF(j))/dr^2+Y(6)*r(j)*cD(j));
        dGdt(j)=De(7)/Ep*((2*cG(j+1)-2*cG(j))/dr^2+Y(7)*r(j)*cD(j));
        dHdt(j)=De(8)/Ep*((2*cH(j+1)-2*cH(j))/dr^2+Y(8)*r(j)*cD(j));
        dIdt(j)=De(9)/Ep*((2*cI(j+1)-2*cI(j))/dr^2+Y(9)*r(j)*cD(j));
    elseif j==N
        
        dAdt(j)=-3/alfa*De(1)/R*((cA(j)-cA(j-1))/dr+Y(1)*r(j)*cD(j))+D*(Feed(1)-cA(N));
        dBdt(j)=-3/alfa*De(2)/R*((cB(j)-cB(j-1))/dr+Y(2)*r(j)*cD(j))+D*(Feed(2)-cB(N));
        dCdt(j)=-3/alfa*De(3)/R*((cC(j)-cC(j-1))/dr+Y(3)*r(j)*cD(j))+D*(Feed(3)-cC(N));
        dDdt(j)=-3/alfa*De(4)/R*((cD(j)-cD(j-1))/dr+r(j)*cD(j))+D*(Feed(4)-cD(N));
        dEdt(j)=-3/alfa*De(5)/R*((cE(j)-cE(j-1))/dr+Y(5)*r(j))*cD(j)+D*(Feed(5)-cE(N));
        dFdt(j)=-3/alfa*De(6)/R*((cF(j)-cF(j-1))/dr+Y(6)*r(j))*cD(j)+D*(Feed(6)-cF(N));
        dGdt(j)=-3/alfa*De(7)/R*((cG(j)-cG(j-1))/dr+Y(7)*r(j))*cD(j)+D*(Feed(7)-cG(N));
        dHdt(j)=-3/alfa*De(8)/R*((cH(j)-cH(j-1))/dr+Y(8)*r(j))*cD(j)+D*(Feed(8)-cH(N));
        dIdt(j)=-3/alfa*De(9)/R*((cI(j)-cI(j-1))/dr+Y(9)*r(j))*cD(j)+D*(Feed(9)-cI(N));
        

    else
        
        dAdt(j)=De(1)/Ep*((cA(j+1)-2*cA(j)+cA(j-1))/dr^2+2/rp(j)*(cA(j)-cA(j-1))/dr+Y(1)*r(j)*cD(j));
        dBdt(j)=De(2)/Ep*((cB(j+1)-2*cB(j)+cB(j-1))/dr^2+2/rp(j)*(cB(j)-cB(j-1))/dr+Y(2)*r(j)*cD(j));
        dCdt(j)=De(3)/Ep*((cC(j+1)-2*cC(j)+cC(j-1))/dr^2+2/rp(j)*(cC(j)-cC(j-1))/dr+Y(3)*r(j)*cD(j));
        dDdt(j)=De(4)/Ep*((cD(j+1)-2*cD(j)+cD(j-1))/dr^2+2/rp(j)*(cD(j)-cD(j-1))/dr+r(j)*cD(j));
        dEdt(j)=De(5)/Ep*((cE(j+1)-2*cE(j)+cE(j-1))/dr^2+2/rp(j)*(cE(j)-cE(j-1))/dr+Y(5)*r(j)*cD(j));
        dFdt(j)=De(6)/Ep*((cF(j+1)-2*cF(j)+cF(j-1))/dr^2+2/rp(j)*(cF(j)-cF(j-1))/dr+Y(6)*r(j)*cD(j));
        dGdt(j)=De(7)/Ep*((cG(j+1)-2*cG(j)+cG(j-1))/dr^2+2/rp(j)*(cG(j)-cG(j-1))/dr+Y(7)*r(j)*cD(j));
        dHdt(j)=De(8)/Ep*((cH(j+1)-2*cH(j)+cH(j-1))/dr^2+2/rp(j)*(cH(j)-cH(j-1))/dr+Y(8)*r(j)*cD(j));
        dIdt(j)=De(9)/Ep*((cI(j+1)-2*cI(j)+cI(j-1))/dr^2+2/rp(j)*(cI(j)-cI(j-1))/dr+Y(9)*r(j)*cD(j));
        
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

dQ=ksis*abs(dH);

dxdt=[dAdt dBdt dCdt dDdt dEdt dFdt dGdt dHdt dIdt dVr dQ]';