clear
close all
load cri1
 Ga1=Ga1_old;

% [Gacrit,Trcrit,Trc,Grc]=Gam_critUall(Tstor,Ga1,Mcrit,ficri);
 [Gacrit,Trcrit,Trc,Grc]=Gam_critScatt(Tstor,Ga1,Mcrit,ficri,fmlst);
 
 map(log10(abs(Gacrit)))