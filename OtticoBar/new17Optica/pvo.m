close all
clear all
%load sap
load sac
clear C D

TinizioF=clock;


mat_mix_vortex

Tfine=clock;

Elt=etime(Tfine,TinizioF)


return

 K1=Kospd;
 K2=Kosmd;
 K3=Koszpd;
 K4=Koszmd;
 save sa K1 K2 K3 K4

clear K1 K2 K3 K4
 K1{11}=Kosp;
 K2{11}=Kosm;
 K3{11}=Koszp;
 K4{11}=Koszm;
 save sa K1 K2 K3 K4

