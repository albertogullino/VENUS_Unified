close all
clear
istrumix=1;
aloc=3.75;
bloc=2.75;
bloc=3.75;
alim=.1;
kcav=28;
kcav0=28;
ifp=4;
dap=.0001;
%dap=.2;
a=30;
adis=30;
igint=1;
nubesu=2;
pimu=[1 2 3];
mbv=[0 1 2];
meun=0;
numodiacc=3;
pasnu=1;
ipolar=1;

cce=.2+i*0.1;
cce=.5+i*1;
%cce=0;
Pshi{1}=4;
Pshi{2}(1)=aloc;
Pshi{3}(1)=bloc;
%Pshi{3}(1)=aloc;
Pshi{4}(1)=dap;
Pshi{5}(1)=cce;

s_sce
