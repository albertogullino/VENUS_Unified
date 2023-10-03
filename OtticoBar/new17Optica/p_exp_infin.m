%load saP
%load saPg
clear 
close all
%load scas
load scainf
%load saPg

pu=[1 21]+19;

Psav=P;

P=Psav(pu,pu);

pausak

M=expm(P);
Ta=M;
lve=length(M)/2;
l1=1:lve;
l2=l1+lve;
Ta11=Ta(l1,l1);
Ta12=Ta(l1,l2);
Ta21=Ta(l2,l1);
Ta22=Ta(l2,l2);
iT=inv(Ta11);

s22=-iT*Ta12;
s21=iT;
s12=Ta22-Ta21*iT*Ta12;
s11=Ta21*iT;

Tscattuu

nv=nitr(1);
fn=(rr/nv);
fr=sqrt(1-(fn*KK).^2);
bea=kcav/fn.*fr*dos;

kzme=kcav/fn;
kt=KK*kcav;
bez=dos*sqrt(kzme^2-kt.^2);


figure            
plot(diag(abs(M-xTras))) 
keyboard

map(log10(abs(M-xTras))) 
title(' differenza matrici trasmissione')
pausak
map(log10(abs(s11))) 
title(' s11 numerico')
pausak
map(log10(abs(s11u))) 
title(' s11 scatt ')

return
, pausak

map(abs(Oo-Tu))            , pausak