function De=DiffusionCoef(T,P,C,i)
Ep=0.703;
Tau=1.17;

if i==7||i==8
        C=[C(7) C(8)];
        D=DG(T,P,C,i);
else
    if i>4&&i<7
        i=i-1;
    end
        C=[C(1) C(2) C(3) C(5) C(6)];
        D=DL(T,P,C,i);
end
De=D*Ep/Tau;

function D=DL(T,P,C,i)
V=1;
MW=[180.156 150.13 80.043 88.11 60.052]; % 8% popol
x=C.*V./MW./sum(C.*V./MW);
Vc=[0.46 0.39 0.16 0.28 0.17];
Pc=[6631.37 6588.38	7355.34 4000 5786.00]/100;
Tc=[1034.02 900.63	677.44 624.00 592.70];
Zc=Pc.*Vc./Tc/8.314;
romol=sum(C.*x)/sum(x)/(sum(MW.*x)/sum(x));
xA=x(i);
x(i)=0;
xB=1-xA;
MW_B=sum(MW.*x)/sum(x);
Pc_B=sum(Pc.*x)/sum(x);
Vc_B=sum(Vc.*x)/sum(x);
Tc_B=sum(Tc.*x)/sum(x);
Zc_B=sum(Zc.*x)/sum(x);

psatunit=["Pa","Pa","mmHg","bar","bar","bar"];
A=[53.16 46.29  10.708 6.10954 	4.68206];
B=[23382 19006  4670 2634.471 1642.54];
C=[0 0 0 -3.471 -39.764];

Tr=0.7;
Ts=Tc.*Tr;
Ps=10.^(A-B./(C+Ts));
Ps(1)=exp(A(1)-B(1)./(C(1)+Ts(1)));
Ps(2)=exp(A(2)-B(2)./(C(2)+Ts(2)));
for j=1:length(Ps)
    switch(psatunit(j))
        case "Pa"
            Ps(j)=Ps(j)*1e-5;
        case "mmHg"
            Ps(j)=Ps(j)/750.061683;
    end
end
    
Pr=Ps./Pc;
omega=-log10(Pr)-1;
omega_B=sum(omega.*x)/sum(x);

sigA=1.866*Vc(i)^(1/3)*Zc(i)^(-6/5);
sigB=1.866*Vc_B^(1/3)*Zc_B^(-6/5);
sigAB=(sigA+sigB)/2;

epsA=65.3*Tc(i)*Zc(i)^(18/5);
epsB=65.3*Tc_B*Zc_B^(18/5);
epsAB=(epsA*epsB)^0.5;
Tx=T/epsAB;
Om=1.06036/Tx^0.1561+0.193*exp(-3.89411*Tx)+1.03587*exp(-1.52996*Tx);

roD0=2.2648*10^-6*T^0.5*(1/MW(i)+1/MW_B)^0.5/sigAB^2/Om;

TrA=T/Tc(i);
TrB=T/Tc_B;
EpA=Tc(i)^(1/6)/(MW(i)^0.5*(0.987*Pc(i))^(2/3));
EpB=Tc_B^(1/6)/(MW_B^0.5*(0.987*Pc_B)^(2/3));

if TrA<1.5
    viskA0=34e-5*TrA^0.94/EpA;
else
    viskA0=17.78e-5*(4.58*TrA-1.67)^(5/8);
end

if TrB<1.5
    viskB0=34e-5*TrB^0.94/EpB;
else
    viskB0=17.78e-5*(4.58*TrB-1.67)^(5/8)/EpB;
end

visk0=(xA*viskA0*MW(i)^0.5+xB*viskB0*MW_B)/(xA*MW(i)^0.5+xB*MW_B^0.5);

romolr=(xA*Vc(i)+xB*Vc_B)*romol;
Eptot=(xA*Tc(i)+xB*Tc_B)^(1/6)/((xA*MW(i)+xB*MW_B)^0.5*(0.987*(xA*Pc(i)+xB*Pc_B))^(2/3));
visk=((0.1023+0.023364*romolr+0.058533*romolr^2-0.040758*romolr^3+0.0093324*romolr^4)^4+10^-4)/Eptot+visk0;

omegaTot=xA*omega(i)+xB*omega_B;
PcTot=xA*Pc(i)+xB*Pc_B;

a=1.07;
b=-0.27-0.38*omegaTot;
c=-0.05+0.1*omegaTot;
PrTot=P/PcTot;

D=a*(visk/visk0)^(b+c*PrTot)*roD0/romol;


function D=DG(T,P,C,i)
MW=[2.002 44.01];
v=[2.31*2 15.9+6.11*2];
MAB=2*(1/MW(1)+1/MW(2))^(-1);
D=0.00143*(T+273.15)^1.75/(P*MAB^0.5*(v(1)^(1/3)+v(2)^(1/3))^2)*1e-4;