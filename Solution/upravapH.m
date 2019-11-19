function n=upravapH(V,cBA,currentpH,targetpH)
K_BA=4.82;
cm=cBA/88.11;

if currentpH>targetpH
    AB='A';
else
    AB='B';
end

switch AB

    case 'A'
        pH1=currentpH;
        pH2=targetpH;
        n=0.5*V*((10^-pH1-10^-pH2)+cm*(10^-pH2/(10^-pH2+10^-K_BA)-10^-pH1/(10^-pH1+10^-K_BA)));
    case 'B'
        pH1=targetpH;
        pH2=currentpH;
        n=V*((10^-pH1-10^-pH2)+cm*(10^-pH2/(10^-pH2+10^-K_BA)-10^-pH1/(10^-pH1+10^-K_BA)));
end
