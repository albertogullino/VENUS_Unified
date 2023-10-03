clear all
close all

lambda=1.4987;
period=.1;
%period=.5;
DC=.5;
Li=.1;
%Li=.001;
iret_BW=0;
itetm=1;
r_in=1;
r_out=3.12;


  ngrating=[4 1];
  par_grat.n1=ngrating(1);
  par_grat.n2=ngrating(2);
  nr1=ngrating(1);
  nr2=ngrating(2);
  par_grat.itetm=itetm;
  par_grat.px=period;
  par_grat.DC=DC;
  NModi=11;
  par_grat.NModi=NModi;
  par_grat.iret_BW=iret_BW;
  par_grat.itetm=itetm;
  par_grat.r_in=r_in;
  par_grat.r_out=r_out;  
  
  phii=30;
  thetai=15;
  Tcampi(lambda,Li,par_grat,phii,thetai);  