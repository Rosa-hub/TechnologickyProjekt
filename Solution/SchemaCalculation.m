function SchemaCalculation

mToExt=2500*17;
Yex1=0.95;
mextOut=2500;
tol=1;

while tol>0.1

  Ferm=Fermentation(mToExt,Yex1,mextOut);
  Ext1=pokus_o_kompletnost_m_TP_dry(Ferm.BA_ext,Ferm.V_ext);
  
  tol=round(abs(Yex1*100-Ext1.Y),1);
  Yex1=Ext1.Y/100; 
  mextOut=Ext1.m_vodna_out;

end
  Ferm=Fermentation(mToExt,Yex1,mextOut);
  Ext1=pokus_o_kompletnost_m_TP_dry(Ferm.BA_ext,Ferm.V_ext);

  disp(Ferm)
  disp(Ext1)
%   asp=AspInteraction