function Tgrating(lambda,thick,par_grat,phii,tetai)

  if ~exist('Nx_in')
   Nx=1;
  else
   Nx=Nx_in;
  end
  neq=0;
  r_in=par_grat.r_in;
  r_out=par_grat.r_out;
  r1=par_grat.n1;
  r2=par_grat.n2;
  period=par_grat.px;
  DC=par_grat.DC;
  NModi=par_grat.NModi;
  d1=period*DC;
  d2=period*(1-DC);
  ite=par_grat.itetm;
  iret_BW=0;
  if isfield(par_grat,'iret_BW')==1
   iret_BW=par_grat.iret_BW;
  end
  Li=thick;
  


  [T11,T22,T12,T21,s11]=orta_skewc(phii,tetai,r_in,r_out,r1,r2,d1,d2,thick,lambda,NModi,0);

