function output=pokus_o_kompletnost_m_TP_dry(c_ferm,V_ferm,mBAinput)
clc
close all
%% pocita pre konst. parametre vsetko
%% vstup MB
% molove hmotnosti
M_BA = 88.11; %kg/kmol
M_h2o = 18.02; %kg/kmol
M_dodecane = 170.34; %kg/kmol
M_IL = 773.27; %kg/kmol 
% vstupne udaje
V_broth =V_ferm; % 10 m3/h % MENIME
c_c0 = c_ferm;   %kmol/m3 % MENIME ferm medium
% c_d0 = 0.0801; %kmol/m3 % MENIME rozpustadlo


m_BA = c_c0*M_BA*V_broth; %kg/h
hust_broth = (0.99704+0.00499421*c_c0-0.000986381*c_c0^2)*1000; % korelacka od martaka, hustota vodneho roztoku BA pri roznych konc kg/m3
m_h2o = V_broth*hust_broth-m_BA; % kg/h % menime
m_vodna_in = m_BA+m_h2o;


m_BA_organika_vstup =mBAinput;% xlsO.xlsDataRead('N43');%xlsread('Connection\MB_dry.xlsx','DM','N43') %kg/h ZADAJ


% organika = IL+DOD+WATER+BA
m_IL = 450; %kg/h              % MENIME ionovka
w_IL0 = 0.6996; %-                  
w_h2o_v_IL0 = 0.0004; %-
w_dodecane0 = 1-w_IL0-w_h2o_v_IL0; %-
m_org0 = m_IL/w_IL0; %kg/h

m_org = m_org0+m_BA_organika_vstup; %kg/h
w_IL = m_IL/m_org; %-                  
w_h2o_v_IL = w_h2o_v_IL0*m_org0/m_org; %-
w_dodecane = w_dodecane0*m_org0/m_org; %-
w_BA_org_vstup = 1-w_IL-w_h2o_v_IL-w_dodecane; %-

m_voda_v_org = w_h2o_v_IL*m_org; %kg/h
m_dodecane = w_dodecane*m_org; %kg/h

w_BA = m_BA/(m_BA+m_h2o); %-
w_h2o = 1-w_BA; %-
n_IL = m_IL/M_IL; %kmol/h
n_dodecane = m_dodecane/M_dodecane; %kmol/h
% hustoty a viskozity 
hustota_h2o_25 = 997; %kg/m3
hustota_BA_25 = 952.8; %kg/m3
hustota_pure_IL_25 = 884.55; %kg/m3
hustota_pure_dodecane_25 = 746.18; %kg/m3
% vodna faza
ro_c = 1/(w_h2o/hustota_h2o_25+w_BA/hustota_BA_25); %kg/m3
vis_c = 0.894e-3;          % Pa.s      % visk. vody
% organicka faza
hustota = [746.18	760.19	784.25	810.99	825.02	839.75	856.22	880.31	885.90	884.55]; %kg/m3
vis = [1.35  1.80	4.90	19.42	26.68	40.18	90.94	244.10	505.01	1058.22]*1e-3; %Pa.s
w_IL_tab = [0	10.10	29.96	50.15	59.94	69.98	79.72	89.81	94.37	99.94]/100; %-
vis_d = interp1(w_IL_tab,vis,w_IL,'spline'); %Pa.s
ro_d = interp1(w_IL_tab,hustota,w_IL,'spline'); %kg/m3

d_ro = abs(ro_d-ro_c);      % kg/m3     % rozdiel hustot
% povrchove napatia - PN 25°C latka - vzduch [mN/m]
PN_IL = 28.32;
PN_dodecane = 24.908;
PN_h2o = 72;
PN_BA = 26.19;
% zmesne PN [mN/m]
% vodna
PN_vodna = PN_BA*w_BA + PN_h2o*w_h2o;
% organicka
PN_org = PN_IL*w_IL + PN_dodecane*w_dodecane + PN_h2o*w_h2o_v_IL ;
% medzifazove napatie [N/m]
MN = (PN_vodna - PN_org)*1e-3;
%% charakteristika kolony
Dc = 0.6;                     % m  1      % priemer kolony
f = 4;                      % s-1 4      % frekvencia
A = 0.02;                   % m         % amplituda
Af = A*f;                   % m.s-1     % rychlost miesania
% konstanta
g = 9.81;                   % m.s-2     % gravitacne zrychlenie
%% rychlosti prudenia faz %m/s
Vc = V_broth/3600*4/pi/Dc^2;         % m/s  % rychlost prudenia kontinalnej fazy     
Vd = m_org/ro_d/3600*4/pi/Dc^2;    % m/s % rychlost prudenia dispergovanej fazy
%% zadrz %-
xd = (3.25e3+7.54e7*Af^3)*Vd^0.81*(Vc+Vd)^0.32*d_ro^-0.98;
%% priemer kvapky %m
C_F = 0.95; % prestup k --> d
C_O = 1.3;
C_P = 0.67;
N = 0.5;
S = 0.5719; % volny prierez platne % m2 % 0.5-0.6 teoreticky 
hc = 0.05; % m
co = 0.7;  % 0.6 kumar alebo 0.7 kathryn orifice coefficient
fi = 2*pi^2/3*((1-S^2)/(hc*co^2*S^2))*Af^3; %disipovana mechanicka energia W/kg
d32 = C_F*S^N*(1/(C_O*(MN/d_ro/g)^0.5)^2+1/(C_P*fi^-0.4*(MN/ro_c)^0.6)^2)^-0.5 ;%m

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
Vs = (-9e-3+2.36e-2*exp(-8.39*Af))*Vd^0.25*MN^0.18*d_ro^0.7;
V = Vd/xd+Vc/(1-xd);
v_terminal = 35.5*(ro_c/1000)^-0.45*(d_ro/1000)^0.58*(vis_c*10)^-0.11*(d32*100)^0.7/100; %m/s

%% D_disp
V_org = m_org/ro_d ;%m3/h
c_d0 = m_BA_organika_vstup/M_BA/V_org; %kmol/m3 % MENIME rozpustadlo
c_IL = n_IL/V_org; %kmol/m3
c_dodecane = n_dodecane/V_org; %kmol/m3
fi_IL_tab = 10.1; %- %25°C
fi_MO = (c_IL*M_IL/hustota_pure_IL_25)*fi_IL_tab+c_dodecane*M_dodecane/hustota_pure_dodecane_25; %-
M_MO = (c_IL*M_IL+c_dodecane*M_dodecane)/(c_IL+c_dodecane); %kg/kmol
vis_MO = vis_d*1e3;          % mPa.s = cP  
V_m_21 = 1388.2; %cm3/mol vypocita excel, komplex (2*BA,1*IL,2*H2O)
T = 25+273.15;       %K
D_disp = 7.4e-8*(fi_MO*M_MO)^0.5*T/(vis_MO*(V_m_21^0.6))*1e-4; %m2/s
%% PL kumar 
Nn = 0.753;
kox_a = Vc^Nn*Af*(1+0.035*fi^0.4*(Vc/Vd)^0.1);
kox_stella = kox_a/a;
dh = 0.0135;
koxa_joseph = 0.43*Af^0.84*dh^-0.21*S^-0.44*hc^-0.41*Vd^0.91;
koxjoseph = koxa_joseph/a;
nastrel_1 = 1e-6;
n_2 = 1/3;
c_2 = 2.44;
k_disp_kumar = fsolve(@PL1,nastrel_1,[opt],d32,Vs,ro_c,ro_d,vis_c,vis_d,fi,g,c_2,n_2,MN,D_disp); %m/s
function f1 = PL1(kd,d32,Vs,ro_c,ro_d,vis_c,vis_d,fi,g,c_2,n_2,MN,D_disp)

Shd = kd*d32/D_disp;
Re = d32*Vs*ro_c/vis_c;
K = vis_d/vis_c;
f1 = 17.7+(3.19e-3*(Re*Shd^(1/3))^1.7)/(1+1.43e-2*(Re*Shd^(1/3))^0.7)*(ro_d/ro_c)^(2/3)*1/(1+K^(2/3))*(1+c_2*(fi/g*(ro_c/g/MN)^(1/4))^n_2)-Shd;
end

nastrel_2 = 1e-7;
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
alfa_c = ((1+1/qc)*exp(Pe_c)-1)^-1; % po porovnani s grafmi to ciselne sedi aj logicky backmixing coefficient between adjacent stages
% ked sa alfa_c = 0, konc su o malinko vyssie coz je OK
%% MB
if exist('n_input','var')==0
    n = 80;     % n celk pocet etazi 
else
    n=n_input;
end
Fc = V_broth/3600; %m3/s
Fd = m_org/(3600*ro_d); %m3/s
S_kolony = pi*Dc^2/4; %m2

c_kont = linspace(c_c0,0,n);
c_disp = linspace(2,c_d0,n);

nastrel = [c_kont' c_disp'];
riesenie = fsolve(@problem_MB,nastrel,[opt],n,c_c0,c_d0,c_IL,k_kont_kumar,k_disp_kumar,a,S_kolony,hc,Fc,Fd,alfa_c);
c_kont = riesenie(:,1);
c_disp = riesenie(:,2);
etaze = 1:n;


vytazok_BA = (c_disp(1)-c_d0)*Fd/(Fc*c_c0)*100 ;% v perc.%
c_kont_out = c_kont(n); %kmol/m3
c_disp_out = c_disp(1);

% 
% figure(1)
% plotyy(etaze,c_kont,etaze,c_disp)

function f = problem_MB(x,n,c_c0,c_d0,c_IL,k_kont_kumar,k_disp_kumar,a,S_kolony,hc,Fc,Fd,alfa_c)
c_c = x(:,1);
c_d = x(:,2);

for i = 1:n
D_1(i) = -642.7879*c_c(i)^3 + 572.6831*c_c(i)^2 -184.8948*c_c(i)+28.8567;
m(i) = D_1(i)*c_IL/0.72;
K_oc_kumar(i) = 1/(1/k_kont_kumar+1/(m(i)*k_disp_kumar)); 

    if i == 1
        f(i,1) = c_c0-(1+alfa_c)*c_c(i)+c_c(i+1)*alfa_c-K_oc_kumar(i).*a*S_kolony*hc/Fc*(c_c(i)-c_d(i)./m(i));
        f(i,2) = c_d(i+1)-c_d(i)+K_oc_kumar(i).*a*S_kolony*hc/Fd*(c_c(i)-c_d(i)./m(i)); % c_d(1) = c_d_out
    elseif i == n
        f(i,1) = (1+alfa_c)*c_c(i-1)-(1+alfa_c)*c_c(i)-K_oc_kumar(i).*a*S_kolony*hc/Fc*(c_c(i)-c_d(i)./m(i)); %c_c(n) = c_c_out
        f(i,2) = c_d0-c_d(i)+K_oc_kumar(i).*a*S_kolony*hc/Fd*(c_c(i)-c_d(i)./m(i));
    else
        f(i,1) = (1+alfa_c)*c_c(i-1)-(1+2*alfa_c)*c_c(i)+alfa_c*c_c(i+1)-K_oc_kumar(i).*a*S_kolony*hc/Fc*(c_c(i)-c_d(i)./m(i));
        f(i,2) = c_d(i+1)-c_d(i)+K_oc_kumar(i).*a*S_kolony*hc/Fd*(c_c(i)-c_d(i)./m(i));
    end

end
end
%% voda v org. faze (IL) na vystupe z extraktora
x =[
0
0
0.005581395
0.00744186
0.020465116
0.044651163
0.106046512
0.195348837
0.31627907
0.524651163
0.641860465
0.744186047];
y = [4.697247706
4.458715596
3.321100917
2.220183486
1.174311927
0.71559633
0.697247706
0.71559633
0.95412844
1.064220183
1.082568807
1.247706422];

%% zadaj:
% Vorg = 0.8060; %m3/h
% c_IL= n_IL/Vorg % bude treba zadat co pouzijeme v extraktore
% c_IL = 0.8022;
CAF = c_c0; % vyplynie z vstupnej konc. BA do fermentora, kt. je v rovnovahe s vystupujucou konc. BA org. fazy 
caf = linspace(0,0.8,100);
for j = 1:length(caf)
if CAF <= 0.045
   cws(j) = -104229*caf(j)^3+10326*caf(j)^2-340.23*caf(j)+4.5896;
else
   cws(j) = 0.7799*caf(j)+0.6386;
end
end
c_IL_exp = 0.64; %kmol/m3
cws_pri_c_IL_exp = interp1(caf,cws,CAF);
cws_pri_c_IL = cws_pri_c_IL_exp*c_IL/c_IL_exp; 
% cws_pri_c_IL  [kmol/m3]to je to co hladame - rovnovazna koncentracia vody
% v rozpustadle pre nasu konc. IL, v rovnovahe s vodnou fazou na vstupe do
% extraktora
% podla mna cws_pri_c_IL = n_vody/Vorg
% cW,S = equilibrium concentration of water in the solvent, kmol·m?3

m_vody = M_h2o*cws_pri_c_IL*V_org-m_voda_v_org;  %kg/h
%% MB vystup
m_BA_PL_z_cont = (c_c0-c_kont_out)*Fc*3600*M_BA ;%kg/h
m_BA_PL_do_disp = (c_disp(1)-c_d0)*Fd*3600*M_BA ;%kg/h
m_organika_out = m_org+m_BA_PL_z_cont+m_vody;
w_BA_org = (m_BA_PL_do_disp+m_BA_organika_vstup)/m_organika_out;
w_dod_org = m_dodecane/m_organika_out;
w_IL_org = m_IL/m_organika_out;
w_h2o_org_out = (m_vody+m_voda_v_org)/m_organika_out;
m_vodna_out = m_vodna_in-m_BA_PL_z_cont-m_vody;
w_BA_vodna = (m_BA-m_BA_PL_z_cont)/m_vodna_out;
w_h2o_vodna = 1-w_BA_vodna;

chyba_MB = m_org+m_vodna_in-m_vodna_out-m_organika_out;

% str1 = sprintf('m_organika_in = %0.1f kg/h \n w_IL = %0.4f \n w_dod = %0.4f \n w_h2o = %0.4f \n w_BA = %0.4f',m_org,w_IL,w_dodecane,w_h2o_v_IL,w_BA_org_vstup);
% disp(str1);
% str2 = sprintf('m_vodna_in = %0.1f kg/h \n w_BA = %0.4f \n w_h2o = %0.4f',m_vodna_in,w_BA,w_h2o);
% disp(str2);
% str3 = sprintf('m_organika_out = %0.1f kg/h \n w_BA = %0.4f \n w_dod = %0.4f \n w_IL = %0.4f \n w_h2o = %0.4f',m_organika_out,w_BA_org,w_dod_org,w_IL_org,w_h2o_org_out);
% disp(str3);
% str4 = sprintf('m_vodna_out = %0.1f kg/h \n w_BA = %0.4f \n w_h2o = %0.4f',m_vodna_out,w_BA_vodna,w_h2o_vodna);
% disp(str4);

xlsO=xlsInteraction;
xlsO.file='MB_dry.xlsm';
xlsO.sheet='EXT';
xlsO=xlsO.xlsOpenConnection;
xlsO.xlsDataWrite('B15',m_vodna_in);
xlsO.xlsDataWrite('B17',w_BA);
xlsO.xlsDataWrite('B19',w_h2o);
xlsO.xlsDataWrite('C15',m_vodna_out);
xlsO.xlsDataWrite('C17',w_BA_vodna);
xlsO.xlsDataWrite('C19',w_h2o_vodna);
xlsO.xlsDataWrite('D15',m_org);
xlsO.xlsDataWrite('D17:D20',[w_BA_org_vstup w_dodecane w_h2o_v_IL w_IL]');
xlsO.xlsDataWrite('E15',m_organika_out);
xlsO.xlsDataWrite('E17:E20',[w_BA_org w_dod_org w_h2o_org_out w_IL_org]');
xlsO.xlsDataWrite('B2:B13',[V_broth c_c0 m_BA_organika_vstup m_IL w_IL0 w_h2o_v_IL0 Dc n vytazok_BA xd xdf a]');
xlsO=xlsO.xlsCloseConnection;
clear xlsO

output.Y=vytazok_BA;
output.m_org_out=m_organika_out;
output.wBA_org=w_BA_org;
output.wdod_org=w_dod_org;
output.wh2o_org=w_h2o_org_out;
output.wIL_org=w_IL_org;
output.m_vodna_in=m_vodna_in;
output.wBA=w_BA;
output.wh2o=w_h2o;
output.m_vodna_out=m_vodna_out;
output.wBA_vodna=w_BA_vodna;
output.wh2o_vodna=w_h2o_vodna;
output.m_org=m_org;
output.wBA_org_vstup=w_BA_org_vstup;
output.wdod=w_dodecane;
output.w_h2o_vIL = w_h2o_v_IL;
output.wIL = w_IL;
output.V_broth=V_broth;
output.c_c0=c_c0;
output.mBA_org_vstup=m_BA_organika_vstup;
output.mIL = m_IL;
output.w_IL0=w_IL0;
output.wh2o_IL0=w_h2o_v_IL0;
output.Dc=Dc;
output.n=n;
output.Hc=n*hc;
output.xd=xd;
output.xdf=xdf;
output.a=a;



caf_I = linspace(0,0.045,100);
for k = 1:length(caf_I)
    cws_I(k) = -104229*caf_I(k)^3+10326*caf_I(k)^2-340.23*caf_I(k)+4.5896;
end 

caf_II = linspace(0.045,0.8,100);
for k = 1:length(caf_II)
    cws_II(k) = 0.7799*caf_II(k)+0.6386;
end
% figure(2)
% plot(caf_I,cws_I,caf_II,cws_II,x,y,'o',CAF,cws_pri_c_IL_exp,'d')
% title('c(ws) = f(c(Af))')
% xlabel('c(Af) [kmol/m3]')
% ylabel('c(ws) [kmol/m3]')
% legend('n = 3','n=1','exp','hladane')
%    
end
