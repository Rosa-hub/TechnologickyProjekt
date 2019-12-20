function SchemaCalculation
clc
clear all
tic
mToExt=2500*17;
Yex1_n=0.95;
mextOut=2500;
mBA_nastrel=0;
tolmain=1;
Y1opt=[];
Y2opt=[];

asp=AspInteraction;
asp=asp.openAspConnection;

while tolmain>0.1
    tol=1;
    while tol>0.1
        Ferm=Fermentation(mToExt,Yex1_n,mextOut);
        Ext1=pokus_o_kompletnost_m_TP_dry(Ferm.BA_ext,Ferm.V_ext,mBA_nastrel);
  
        tol=round(abs(Yex1_n*100-Ext1.Y),1);
        Yex1_n=Ext1.Y/100; 
        mextOut=Ext1.m_vodna_out;

    end



%Input dát do aspenu
asp=asp.clearParams;
asp.propType='TOTFLOW';
asp.streamName='F0';
asp.aspSetScalar(Ext1.m_org_out);
asp.propType='FLOW';
asp.componentName='IL';
asp.aspSetScalar(Ext1.wIL_org);
asp.componentName='BA';
asp.aspSetScalar(Ext1.wBA_org);
asp.componentName='H2O';
asp.aspSetScalar(Ext1.wh2o_org);
asp.componentName='DOD';
asp.aspSetScalar(Ext1.wdod_org);
asp.aspRunSimulation;


%Export dát z aspenu
%Dekanter
asp=asp.clearParams;
asp.propType='MASSFLMX';
asp.streamName='D1C';
mF_dekanter=asp.aspGetScalar;
asp.propType='MASSFRAC';
asp.componentName='BA';
wBA_dekanter=asp.aspGetScalar;
asp.componentName='DOD';
wDOD_dekanter=asp.aspGetScalar;

%SPD_MIX
asp=asp.clearParams;
asp.streamName='D2C';
asp.propType='MASSFLOW';
asp.componentName='BA';
mBA_to_mix=asp.aspGetScalar;
asp.componentName='H2O';
mH2O_to_mix=asp.aspGetScalar;
asp.componentName='DOD';
mDOD_to_mix=asp.aspGetScalar;

%SPD_B_to_EXT
asp=asp.clearParams;
asp.streamName='CIL';
asp.propType='MASSFLOW';
asp.componentName='BA';
mBA_to_e1=asp.aspGetScalar;
asp.componentName='IL';
mIL_to_e1=asp.aspGetScalar;
asp.componentName='DOD';
mDOD_to_e1=asp.aspGetScalar;

YSPD=(1-mBA_to_e1/(Ext1.m_org_out*Ext1.wBA_org))*100;

xlsO=xlsInteraction;
xlsO.file='MB_dry.xlsm';
xlsO.sheet='DM';
xlsO=xlsO.xlsOpenConnection;
xlsO.xlsDataWrite('C2:C4',[mF_dekanter wBA_dekanter wDOD_dekanter]');
xlsO.xlsDataWrite('C18:C23',[mBA_to_mix mH2O_to_mix mDOD_to_mix mDOD_to_e1 mIL_to_e1 mBA_to_e1]');
xlsO=xlsO.xlsCloseConnection;

Ext2=pokus_o_kompletnost_m_TP_2;

xlsO=xlsInteraction;
xlsO.file='MB_dry.xlsm';
xlsO.sheet='DM';
xlsO=xlsO.xlsOpenConnection;
BA_e1=xlsO.xlsDataRead('N43');
xlsO=xlsO.xlsCloseConnection;


% disp(Ferm)
% disp(Ext1)
  
  tolmain=abs(round(BA_e1,1)-round(Ext1.mBA_org_vstup,1));
  mBA_nastrel = round(BA_e1,2);
end

AD=adsorpcia(Ext1.m_vodna_out/1000);
% disp(sprintf('Cas simulacie: %1.0f:%1.0f',round(toc/60,0),round(mod(toc,60),0)))
%     x=xlsInteraction;
%     x.file='Economy_test.xlsm';
%     x=x.xlsOpenConnection;
%     x.xlsObject.Visible=1;
%     x.xlsRunMacro('copyTemplate');  
% 
asp.aspCloseAspConnection;
% x=x.xlsCloseConnection;

disp(Ext2)
disp(Ext1)
disp(Ferm)
disp(AD)


% function Ext1=extrakcia1(BA,Vext,Y)
% 
% if isempty(Y)==1
%     Ext1=pokus_o_kompletnost_m_TP_dry(BA,Vext);
%     
% else
%     n=fsolve(@solExt1,80,[],BA,Vext,Y);
%     Ext1=pokus_o_kompletnost_m_TP_dry(BA,Vext,ceil(n));
% end
% 
% function Ext2=extrakcia2(Y)
% 
% if isempty(Y)==1
%     Ext2=pokus_o_kompletnost_m_TP_dry();
%     
% else
%     n=fsolve(@solExt1,80,[],BA,Vext,Y);
%     Ext2=pokus_o_kompletnost_m_TP_2(ceil(n));
% end
% 
% function f=solExt1(n,BA,Vext,Y)
% n=ceil(n);
% Ext1=pokus_o_kompletnost_m_TP_dry(BA,Vext,n);
% f=round(Ext1.Y/100,1)-Y;
% 
% function f=solExt2(n)
% n=ceil(n);
% Ext1=pokus_o_kompletnost_m_TP_dry(BA,Vext,n);
% f=round(Ext1.Y/100,1)-Y;
