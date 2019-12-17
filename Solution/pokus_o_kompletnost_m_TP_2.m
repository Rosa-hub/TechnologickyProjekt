function output =pokus_o_kompletnost_m_TP_2
clc
close all
%% pocita pre konst. parametre vsetko
%% vstup MB
% molove hmotnosti
M_BA = 88.11; %kg/kmol
M_h2o = 18.02; %kg/kmol
M_dodecane = 170.34; %kg/kmol
M_IL = 773.27; %kg/kmol 
%% dekanter

% vstupne udaje
% V_sol = xlsread('MB_dry.xlsx','DM','M28') %  m3/h % menime
% c_d0 =xlsread('MB_dry.xlsx','DM','L28') %kmol/m3
xlsO=xlsInteraction;
xlsO.file='MB_dry.xlsm';
xlsO.sheet = 'DM';
xlsO=xlsO.xlsOpenConnection;
delSol=xlsO.xlsDataRead('F8');

while abs(delSol)>1e-6
    xlsO.xlsRunMacro('DcoefMacro');
    delSol=xlsO.xlsDataRead('F8');
end

V_sol=xlsO.xlsDataRead('M28');
c_d0=xlsO.xlsDataRead('L28');
z1=xlsO.xlsDataRead('I13');
xlsO=xlsO.xlsCloseConnection;

% V_sol = xlsread('MB_wet.xlsx','DM','M28') %  m3/h % menime
% c_d0 =xlsread('MB_wet.xlsx','DM','L28') %kmol/m3
c_c0 = 0.0;   %kmol/m3 % nemenne
hust_sol = 19.042*c_d0+753.6;% hustota org. roztoku BA  odhad %pri roznych konc kg/m3
m_BA = c_d0*M_BA*V_sol; %kg/h
hust_solv = 746.18;%
m_dod = V_sol*hust_solv; %kg/h
m_sol = m_dod+m_BA; % kg/h % menime
w_dod = m_dod/m_sol; % -
w_BA = 1-w_dod; %-
n_dodecane = m_dod/M_dodecane; %kmol/h
m_water = 220; %kg/h              % menime
n_water = m_water/M_h2o; %kmol/h
w_h2o = 1; % -
% hustoty a viskozity 
hustota_h2o_25 = 997; %kg/m3
hustota_BA_25 = 952.8; %kg/m3
hustota_pure_IL_25 = 884.55; %kg/m3
hustota_pure_dodecane_25 = 746.18; %kg/m3
% vodna faza
V_vod = m_water/hustota_h2o_25;%m3/h
vis_c = 0.894e-3;          % Pa.s      % visk. vody
% organicka faza
hustota = [746.18	757.87	803.03	829.13	841.76	858.42	874.56	895.99	903.38	914.40]; %kg/m3
vis = [1.35 1.88 7.67 14.69	21.73 33.14	52.30 90.26	121.27 156.74]*1e-3; %Pa.s
w_IL_tab = [0.00	10.03	29.19	46.91	55.02	63.29	70.93	79.45	82.33	84.65]/100; %-
vis_d = 1.34e-3; % cisty dodekan %interp1(w_IL_tab,vis,w_IL,'spline'); %Pa.s
ro_d = hust_sol; %interp1(w_IL_tab,hustota,w_IL,'spline'); %kg/m3
ro_c = hustota_h2o_25;
d_ro = abs(ro_d-ro_c);      % kg/m3     % rozdiel hustot
% povrchove napatia - PN 25°C latka - vzduch [mN/m]
PN_dodecane = 24.908;
PN_h2o = 72;
PN_BA = 26.19;
% zmesne PN [mN/m]
% vodna
PN_vodna = PN_h2o*w_h2o;
% organicka
PN_org = PN_dodecane*w_dod + PN_BA*w_BA ;
% medzifazove napatie [N/m]
MN = (PN_vodna - PN_org)*1e-3;
%% charakteristika kolony
Dc = 0.2;                     % m  1      % priemer kolony
f = 4;                      % s-1 4      % frekvencia
A = 0.02;                   % m         % amplituda
Af = A*f;                   % m.s-1     % rychlost miesania
% konstanta
g = 9.81;                   % m.s-2     % gravitacne zrychlenie
%% rychlosti prudenia faz %m/s
Vc = V_vod/3600*4/pi/Dc^2;         % m/s  % rychlost prudenia kontinalnej fazy     
Vd = V_sol/3600*4/pi/Dc^2;    % m/s % rychlost prudenia dispergovanej fazy
%% zadrz %-
xd = (2.14e3+1.65e7*Af^3)*Vd^0.81*(Vc+Vd)^0.32*d_ro^-0.98;
%% priemer kvapky %m
C_F = 1.48; % prestup d --> k
C_O = 1.3;
C_P = 0.67;
N = 0.5;
S = 0.5719; % volny prierez platne % m2 % 0.5-0.6 teoreticky 
hc = 0.05; % m
co = 0.7;  % 0.6 kumar alebo 0.7 kathryn orifice coefficient
fi = 2*pi^2/3*((1-S^2)/(hc*co^2*S^2))*Af^3; %disipovana mechanicka energia W/kg
d32 = C_F*S^N*(1/(C_O*(MN/d_ro/g)^0.5)^2+1/(C_P*fi^-0.4*(MN/ro_c)^0.6)^2)^-0.5; %m

%% zahltenie %- Vc>0
nastrel_zahltenie = 0.3;
L = Vd/Vc;
beta = 24*vis_c/(0.53*d32*ro_c);
gama = 4*d32*g*d_ro/(1.59*ro_c);
opt=optimset('Display','off');
xdf = fsolve(@zahltenie,nastrel_zahltenie,[opt],L,beta,gama);
%Vcf m/s
Vcf = ((-beta+(beta^2+4*gama*(1-xdf)/(1+4.56*xdf^0.73))^0.5)*xdf*(1-xdf))/(2*(xdf+L*(1-xdf)));

% str = sprintf('xdf = %1.4f     Vcf = %1.4f    ',xdf,Vcf)
function f1 = zahltenie(xdf,L,beta,gama)
delta = (beta^2+4*gama*(1-xdf)/(1+4.56*xdf^0.73))^0.5;
f1 = (xdf+L*(1-xdf))*((delta-beta)*(1-2*xdf)-2*gama*xdf*(1-xdf)/(delta*(1+4.56*xdf^0.73)^2)*(1+4.56*xdf^0.73+3.33*xdf^-0.27*(1-xdf)))+(delta-beta)*xdf*(1-xdf)*(L-1);
end

if xd>=xdf
    error
end
%% specificka stycna plocha m-1
a = 6*xd/d32; %m-1
%% rychlost padania kvapky m/s
Vs = (-9.82e-3+3.07e-2*exp(-8.39*Af))*Vd^0.25*MN^0.18*d_ro^0.7;
V = Vd/xd+Vc/(1-xd);
v_terminal = 35.5*(ro_c/1000)^-0.45*(d_ro/1000)^0.58*(vis_c*10)^-0.11*(d32*100)^0.7/100; %m/s

%% D_disp
V_org = m_sol/hust_sol; %m3/h
c_dodecane = n_dodecane/V_org; %kmol/m3
fi_dod_tab = 1; %- %25°C neasociovane rozpustadlo markos prednaska
fi_MO = fi_dod_tab; %-
M_MO = M_dodecane; %kg/kmol
vis_MO = vis_d*1e3;          % mPa.s = cP  
V_m_21 = 104.5; %cm3/mol molovy objem kys. maslovej z excelu
T = 25+273.15;       %K
D_disp = 7.4e-8*(fi_MO*M_MO)^0.5*T/(vis_MO*(V_m_21^0.6))*1e-4; %m2/s
%% PL kumar 
Nn = 0.753;
kox_a = Vc^Nn*Af*(1+0.035*fi^0.4*(Vc/Vd)^0.1);
kox_stella = kox_a/a;
dh = 0.0135;
koxa_joseph = 0.43*Af^0.84*dh^-0.21*S^-0.44*hc^-0.41*Vd^0.91;
koxjoseph = koxa_joseph/a;
nastrel_1 = 0.5e-5;
n_2 = 1/3;
c_2 = 2.44;
k_disp_kumar = fsolve(@PL1,nastrel_1,[opt],d32,Vs,ro_c,ro_d,vis_c,vis_d,fi,g,c_2,n_2,MN,D_disp); %m/s
function f1 = PL1(kd,d32,Vs,ro_c,ro_d,vis_c,vis_d,fi,g,c_2,n_2,MN,D_disp)

Shd = kd*d32/D_disp;
Re = d32*Vs*ro_c/vis_c;
K = vis_d/vis_c;
f1 = 17.7+(3.19e-3*(Re*Shd^(1/3))^1.7)/(1+1.43e-2*(Re*Shd^(1/3))^0.7)*(ro_d/ro_c)^(2/3)*1/(1+K^(2/3))*(1+c_2*(fi/g*(ro_c/g/MN)^(1/4))^n_2)-Shd;
end

nastrel_2 = 1e-5;
n_1 = 1/3;
c_1 = 2.44;
k_kont_kumar = fsolve(@PL2,nastrel_2,[opt],d32,Vs,ro_c,vis_c,vis_d,fi,g,c_1,n_1,MN,xd); %m/s
function f2 = PL2(kc,d32,Vs,ro_c,vis_c,vis_d,fi,g,c_1,n_1,MN,xd)
D_kont = 1.0427e-9; %m2/s
Shc = kc*d32/D_kont;
Scc = vis_c/(ro_c*D_kont);
Re = d32*Vs*ro_c/vis_c;
Pec = d32*Vs/D_kont;
Shc_rigid = 2.43+0.775*Re^0.5*Scc^(1/3)+0.0103*Re*Scc^(1/3);
c1 = 50;
Shc_inf = c1+2/(sqrt(pi))*Pec^0.5;
K = vis_d/vis_c;
f2 = 5.26e-2*Re^(1/3+6.59e-2*Re^(1/4))*Scc^(1/3)*(Vs*vis_c/MN)^(1/3)*1/(1+K^1.1)*(1+c_1*(fi/g*(ro_c/g/MN)^(1/4))^n_1)-(Shc/(1-xd)-Shc_rigid)/(Shc_inf-Shc/(1-xd));
end

% m = 14.6019*c_IL/0.72; % - z distribucny_koef.m z clanku martak 2008, precitany graf
% K_od_kumar = 1/(m/k_kont_kumar + 1/k_disp_kumar); %celkovy koef PL v disp faze % ms-1
% K_oc_kumar = 1/(1/k_kont_kumar+1/(m*k_disp_kumar))  %celkovy koef PL v kont faze % ms-1

%% axialna disperzia 
dh = 0.0135; %m %priemer perforacie
M = 2.28*(2*A/pi); % vyska well mixed region %m
Ec_kumar = 0.0688*2*Af*dh/S^1.5+0.082*d32*(4*g*d_ro*(v_terminal-Vc)*xd*(hc-M)/(ro_c*(1-xd)))^1/3; %m2/s
Pe_c = Vc*(hc-M)/Ec_kumar; % -
R = Vc/Vd; %-
FI = asin(Vd*(R-1)/(pi*Af));
% qd nie je 0 lebo emulsion operation so qc sa meni z .. na ..
%qc0 = ((R-1)/(R*pi))*(FI+cot(FI)-(pi/2))-((1+0)/R) %if qd = 0, mixer settler
qc = (1/pi)*((FI+(1-xd)*(R-1)/R)*cot(FI)-(pi/2)); % mozno jedna zatvorka nesedi ale po porovnani s qc0 cca rovnake
alfa_c = 0;%((1+1/qc)*exp(Pe_c)-1)^-1; % po porovnani s grafmi to ciselne sedi aj logicky backmixing coefficient between adjacent stages
% ked sa alfa_c = 0, konc su o malinko vyssie coz je OK
%% MB
if exist('n_input','var')==0
n = 120;     % n celk pocet etazi 

else 
    n=n_input;
end
Fc = V_vod/3600; %m3/s
Fd = V_sol/3600; %m3/s
S_kolony = pi*Dc^2/4; %m2

c_kont = linspace(c_c0,2,n);
c_disp = linspace(0,c_d0,n);

nastrel = [c_kont' c_disp'];
riesenie = fsolve(@problem_MB,nastrel,[opt],n,c_c0,c_d0,k_kont_kumar,k_disp_kumar,a,S_kolony,hc,Fc,Fd,alfa_c);
c_kont = riesenie(:,1);
c_disp = riesenie(:,2);
etaze = 1:n;

c_kont_out = c_kont(n); %kmol/m3
c_disp_out = c_disp(1);%kmol/m3
vytazok_BA = (-c_disp(1)+c_d0)*Fd/(Fd*c_d0)*100;% v perc.%

m_BA_PL_do_cont = (-c_c0+c_kont_out)*Fc*3600*M_BA; %kg/h
m_BA_PL_z_disp = (-c_disp(1)+c_d0)*Fd*3600*M_BA; %kg/h

m_vodna_out = m_BA_PL_z_disp+m_water;
w_BA_vodna = m_BA_PL_z_disp/(m_BA_PL_z_disp+m_water); %-
w_h2o_vodna = 1-w_BA_vodna; %-

m_organika_out = m_dod+m_BA - m_BA_PL_z_disp;
w_BA_org = (m_BA - m_BA_PL_z_disp)/(m_dod+m_BA - m_BA_PL_z_disp);
w_dod_org = 1-w_BA_org;

chyba_MB = m_sol+m_water-m_vodna_out-m_organika_out;

% str1 = sprintf('m_organika_in = %0.1f kg/h \n w_BA = %0.4f \n w_h2o = %0.4f',m_sol,w_BA,w_dod);
% disp(str1);
% str2 = sprintf('m_vodna_in = %0.1f kg/h \n w_BA = 0 \n w_h2o = 1',m_water);
% disp(str2);
% str3 = sprintf('m_organika_out = %0.1f kg/h \n w_BA = %0.4f \n w_h2o = %0.4f',m_organika_out,w_BA_org,w_dod_org);
% disp(str3);
% str4 = sprintf('m_vodna_out = %0.1f kg/h \n w_BA = %0.4f \n w_h2o = %0.4f',m_vodna_out,w_BA_vodna,w_h2o_vodna);
% disp(str4);
% 
% xlswrite('MB_dry.xlsx',m_water,'EXT','J15');
% xlswrite('MB_dry.xlsx',m_vodna_out,'EXT','K15');
% xlswrite('MB_dry.xlsx',w_BA_vodna,'EXT','K17');
% xlswrite('MB_dry.xlsx',w_h2o_vodna,'EXT','K19');
% xlswrite('MB_dry.xlsx',m_sol,'EXT','L15');
% xlswrite('MB_dry.xlsx',w_BA,'EXT','L17');
% xlswrite('MB_dry.xlsx',w_dod,'EXT','L18');
% xlswrite('MB_dry.xlsx',m_organika_out,'EXT','M15');
% xlswrite('MB_dry.xlsx',w_dod_org,'EXT','M18');
% xlswrite('MB_dry.xlsx',w_BA_org,'EXT','M17');
% xlswrite('MB_dry.xlsx',V_sol,'EXT','J2');
% xlswrite('MB_dry.xlsx',c_d0,'EXT','J3');
% xlswrite('MB_dry.xlsx',m_water,'EXT','J4');
% xlswrite('MB_dry.xlsx',Dc,'EXT','J5');
% xlswrite('MB_dry.xlsx',n,'EXT','J6');
% xlswrite('MB_dry.xlsx',vytazok_BA,'EXT','J7');
% xlswrite('MB_dry.xlsx',xd,'EXT','J8');
% xlswrite('MB_dry.xlsx',xdf,'EXT','J9');
% xlswrite('MB_dry.xlsx',a,'EXT','J10');

xlsO.file='MB_dry.xlsm';
xlsO.sheet='EXT';
xlsO=xlsO.xlsOpenConnection;
xlsO.xlsDataWrite('J15',m_water);
xlsO.xlsDataWrite('K15',m_vodna_out);
xlsO.xlsDataWrite('K17',w_BA_vodna);
xlsO.xlsDataWrite('K19',w_h2o_vodna);
xlsO.xlsDataWrite('L15',m_sol);
xlsO.xlsDataWrite('L17:L18',[w_BA w_dod]');
xlsO.xlsDataWrite('L15',m_organika_out);
xlsO.xlsDataWrite('M17:M18',[w_BA_org w_dod_org]');
xlsO.xlsDataWrite('J2:J10',[V_sol c_d0 m_water Dc n vytazok_BA xd xdf a]');
z2=xlsO.xlsDataRead('M28');
xlsO=xlsO.xlsCloseConnection;

output.m_water=m_water;
output.m_vodna_out=m_vodna_out;
output.w_BA_vodna=w_BA_vodna;
output.w_h2o_vodna=w_h2o_vodna;
output.m_sol=m_sol;
output.w_BA=w_BA;
output.w_dod=w_dod;
output.m_organika_out=m_organika_out;
output.w_BA_org=w_BA_org;
output.w_dod_org=w_dod_org;
output.V_sol = V_sol;
output.c_d0 = c_d0;
output.Dc = Dc;
output.n=n;
output.Hc=n*hc;
output.Y=vytazok_BA;
output.xd=xd;
output.xdf=xdf;
output.a=a;
output.z1=z1;
output.z2=z2;

% xlswrite('MB_wet.xlsx',m_water,'EXT','J15');
% xlswrite('MB_wet.xlsx',m_vodna_out,'EXT','K15');
% xlswrite('MB_wet.xlsx',w_BA_vodna,'EXT','K17');
% xlswrite('MB_wet.xlsx',w_h2o_vodna,'EXT','K19');
% xlswrite('MB_wet.xlsx',m_sol,'EXT','L15');
% xlswrite('MB_wet.xlsx',w_BA,'EXT','L17');
% xlswrite('MB_wet.xlsx',w_dod,'EXT','L18');
% xlswrite('MB_wet.xlsx',m_organika_out,'EXT','M15');
% xlswrite('MB_wet.xlsx',w_dod_org,'EXT','M18');
% xlswrite('MB_wet.xlsx',w_BA_org,'EXT','M17');
% xlswrite('MB_wet.xlsx',V_sol,'EXT','J2');
% xlswrite('MB_wet.xlsx',c_d0,'EXT','J3');
% xlswrite('MB_wet.xlsx',m_water,'EXT','J4');
% xlswrite('MB_wet.xlsx',Dc,'EXT','J5');
% xlswrite('MB_wet.xlsx',n,'EXT','J6');
% xlswrite('MB_wet.xlsx',vytazok_BA,'EXT','J7');
% xlswrite('MB_wet.xlsx',xd,'EXT','J8');
% xlswrite('MB_wet.xlsx',xdf,'EXT','J9');
% xlswrite('MB_wet.xlsx',a,'EXT','J10');

% figure(1)
% plotyy(etaze,c_kont,etaze,c_disp)

function f = problem_MB(x,n,c_c0,c_d0,k_kont_kumar,k_disp_kumar,a,S_kolony,hc,Fc,Fd,alfa_c)
c_c = x(:,1);
c_d = x(:,2);

for i = 1:n
D_1(i) = 0.8565*c_c(i) + 0.1609;
m(i) = D_1(i);
K_od_kumar(i) = 1/(m(i)/k_kont_kumar+1/k_disp_kumar); 

    if i == 1
        f(i,1) = c_c0-(1+alfa_c)*c_c(i)+c_c(i+1)*alfa_c-K_od_kumar(i).*a*S_kolony*hc/Fc*(c_c(i)-c_d(i)./m(i));
        f(i,2) = c_d(i+1)-c_d(i)+K_od_kumar(i).*a*S_kolony*hc/Fd*(c_c(i)-c_d(i)./m(i)); % c_d(1) = c_d_out
    elseif i == n
        f(i,1) = (1+alfa_c)*c_c(i-1)-(1+alfa_c)*c_c(i)-K_od_kumar(i).*a*S_kolony*hc/Fc*(c_c(i)-c_d(i)./m(i)); %c_c(n) = c_c_out
        f(i,2) = c_d0-c_d(i)+K_od_kumar(i).*a*S_kolony*hc/Fd*(c_c(i)-c_d(i)./m(i));
    else
        f(i,1) = (1+alfa_c)*c_c(i-1)-(1+2*alfa_c)*c_c(i)+alfa_c*c_c(i+1)-K_od_kumar(i).*a*S_kolony*hc/Fc*(c_c(i)-c_d(i)./m(i));
        f(i,2) = c_d(i+1)-c_d(i)+K_od_kumar(i).*a*S_kolony*hc/Fd*(c_c(i)-c_d(i)./m(i));
    end

end
end

   
end
