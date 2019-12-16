function voda_v_IL
clc
% close all
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
Vorg = 0.8060; %m3/h
n_IL = ; %kmol/h
c_IL= n_IL/Vorg % bude treba zadat co pouzijeme v extraktore
c_IL = 0.8022;
CAF = 0.2; % vyplynie z vstupnej konc. BA do fermentora, kt. je v rovnovahe s vystupujucou konc. BA org. fazy 


caf = linspace(0,1.5,100);
for i = 1:length(caf)
if CAF <= 0.045
   cws(i) = -104229*caf(i)^3+10326*caf(i)^2-340.23*caf(i)+4.5896;
else
   cws(i) = 0.7799*caf(i)+0.6386;
end
end
c_IL_exp = 0.64; %kmol/m3
cws_pri_c_IL_exp = interp1(caf,cws,CAF)
cws_pri_c_IL = cws_pri_c_IL_exp*c_IL/c_IL_exp;
% cws_pri_c_IL  [kmol/m3]to je to co hladame - rovnovazna koncentracia vody
% v rozpustadle pre nasu konc. IL, v rovnovahe s vodnou fazou na vstupe do
% extraktora
% podla mna cws_pri_c_IL = n_vody/Vorg
% cW,S = equilibrium concentration of water in the solvent, kmol·m?3

m_vody = 18.02*cws_pri_c_IL*Vorg %kg/m3



caf_I = linspace(0,0.045,100);
for i = 1:length(caf_I)
    cws_I(i) = -104229*caf_I(i)^3+10326*caf_I(i)^2-340.23*caf_I(i)+4.5896;
end 

caf_II = linspace(0.045,0.8,100);
for i = 1:length(caf_II)
    cws_II(i) = 0.7799*caf_II(i)+0.6386;
end

plot(caf_I,cws_I,caf_II,cws_II,x,y,'o',CAF,cws_pri_c_IL_exp,'d')
title('c(ws) = f(c(Af))')
xlabel('c(Af) [kmol/m3]')
ylabel('c(ws) [kmol/m3]')
legend('n = 3','n=1','exp','hladane')
