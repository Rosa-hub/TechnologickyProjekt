function SchemaCalculation
clc
clear all
tic

Y1opt=[0.90 0.95 0.99];
Y2opt=[0.90 0.95 0.99];

asp=AspInteraction;
asp=asp.openAspConnection;
for i=1:length(Y1opt)
    
    for j=1:length(Y2opt)
       
        mToExt=2500*17;
        Yex1_n=0.95;
        mextOut=2500;
        mBA_nastrel=0;
        tolmain=1;
        
        while tolmain>0.1
            tol=1;
            while tol>0.1
                Ferm=Fermentation(mToExt,Yex1_n,mextOut);
                Ext1=extrakcia1(Ferm.BA_ext,Ferm.V_ext,mBA_nastrel,Y1opt(i));

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
            Ext2=extrakcia2(Y2opt(j));
            
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
            %   tolmain=0;


        end
        
            Y1(i)=Ext1.Y;
            Y2(j)=Ext2.Y;
            Yspd(i,j)=YSPD;
            z1e(i,j)=Ext2.z1;
            z2e(i,j)=Ext2.z2;
        
    end
end
figure(1)
surf(Y1,Y2,z1e)
xlabel 'Y1 [%]'
ylabel 'Y2 [%]'
zlabel 'z1 [%]'

figure(2)
surf(Y1,Y2,z2e)
xlabel 'Y1 [%]'
ylabel 'Y2 [%]'
zlabel 'z2 [%]'


x=xlsInteraction;
x.file='ExportData.xlsx';
x=x.xlsOpenConnection;
x.xlsObject.Visible=1;
x.xlsDataWrite(x.getRange("A1",Y1'),Y1')
x.xlsDataWrite(x.getRange("B1",Y2'),Y2')
x.xlsDataWrite(x.getRange("D1",z1e),z1e)
x.xlsDataWrite(x.getRange(strcat("D",string(length(Y1)+2)),z2e),z2e)
x=x.xlsCloseConnection;

disp(sprintf('Cas simulacie: %1.0f:%1.0f',round(toc/60,0),round(mod(toc,60),0)))
%     x=xlsInteraction;
%     x.file='Economy_test.xlsm';
%     x=x.xlsOpenConnection;
%     x.xlsObject.Visible=1;
%     x.xlsRunMacro('copyTemplate');  

asp.aspCloseAspConnection;
% x=x.xlsCloseConnection;

disp(Ext2)
% disp(YSPD)
disp(Ext1)
disp(Ferm)



function Ext1=extrakcia1(BA,Vext,mBA_nastrel,Y)

if isempty(Y)==1
    
    Ext1=Ext1_opt(BA,Vext,mBA_nastrel);
    
else
    mSOL=fsolve(@solExt1,450,[],BA,Vext,Y,mBA_nastrel);
    Ext1=Ext1_opt(BA,Vext,mBA_nastrel,mSOL);
end

function Ext2=extrakcia2(Y)

if isempty(Y)==1
    Ext2=Ext2_opt;
    
else
    mSOL=fsolve(@solExt2,220,[],Y);
    Ext2=Ext2_opt(mSOL);
end

function f=solExt1(mSOL,BA,Vext,Y,mBA_in)

Ext1=Ext1_opt(BA,Vext,mBA_in,mSOL);
f=(Ext1.Y/100-Y);

function f=solExt2(mSOL,Y)

Ext2=Ext2_opt(mSOL);
f=(Ext2.Y/100-Y);

