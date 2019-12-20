function output = adsorpcia(V_c)
clc 
close all

%% data
rozpustnost = 9.1e-3; % kg/m3 rozpustnost IL vo vode, od martaka
%% experiment
S_kolony = 2e-4; % m2 prierez kolony
m_ADS = 3.77e-3; %kg hmotnost naplne amberlit XAD 1180N
H_naplne = 15.5e-2; % vyska naplne m
pomer = S_kolony/H_naplne; % m
V_naplne = S_kolony*H_naplne; %m3
ro_s = m_ADS/V_naplne; % sypna hustota adsorbenta kg/m3
V_nastreku = 3.8e-6/60; %m3/s objemovy prietok nastreku
v = 1.88e-2/60; %mimovrstvova rychlost m/s zo spravy %nizsia ako vypocitan, pouzivaj vypocitanu
BV = 169.4; % vydrzala napln bez prienikovej krivky
V_celk = BV*V_naplne; %m3  % 5250cm3
w = V_nastreku/S_kolony; %vypoc. mimovrstvova rychlost
t_ala_prienik = V_celk/V_nastreku;  %s max cas exp. a nedoslo k prieniku; na vystupe bez IL (uvedene 22h)
t_p_hod = t_ala_prienik/3600;
mIL = V_nastreku*rozpustnost*t_ala_prienik; %kg hm. IL zachytenej na adsorbente po 23h

%% scale up zadaj V_c a t_ADS
% V_c = 2.65; % m3/h prietok kontinualnej fazy na vystupe z extraktora 
t_ADS = 30*24*3600; % s 

% zvyšovaním prietoku sa zväèšuje priemer kolóny 
% zvyšovaním èasu ADs sa zväèšuje výška kolóny 

V_nastreku_novy = V_c/3600; %m3/s do adsorbera
S_nove = V_nastreku_novy/w; %m2 zostava konst. vyp. medzivrstvova rychlost z exp.
d_nove = sqrt(S_nove*4/pi); %m
mIL_nove = V_nastreku_novy*rozpustnost*t_ADS; %kg
m_ADSnove = mIL_nove*m_ADS/mIL; %kg
V_naplne_nove = m_ADSnove/ro_s; %m3
H_naplne_nove = V_naplne_nove/S_nove; % m

str = sprintf('Pri prietoku kont. f. %0.1f m3/h \na èase adsorpcie %0.1f h \nby bol priemer adsorbera %0.1f \na výška náplne %1.f m. \nNaadsorbuje sa %0.1f kg IL.',V_c,t_ADS/3600,d_nove,H_naplne_nove,mIL_nove);
% disp(str);

%% desorpcia
prietok_desorbenta = S_nove*w; %m3/s
V_des = 12*V_naplne_nove; %m3
t_des = V_des/prietok_desorbenta/3600; %h

output.D=round(d_nove,2);
output.H=round(H_naplne_nove/0.9,1);
output.mIL=mIL_nove;
output.Vdes=prietok_desorbenta*3600;
output.tdes=t_des;


