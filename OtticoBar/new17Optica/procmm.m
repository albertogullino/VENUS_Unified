close all
clear all

load 

ifp=-10
ipcam=1;

unpack_0



     lam0=lambda0*1e6;
     fish=find(shavet==6)
   
     if length(fish)>0  
       ngrating=nv(fish,:);
       par_grat.r_in=nv(fish-1,1);
       par_grat.r_out=nv(fish+1,1);
       par_grat.n1=ngrating(1);
       par_grat.n2=ngrating(2);
       par=radii.array{fish};
       n1g=ngrating(1);
       n2g=ngrating(2);
       period=par{5};
       t1=par{6};
       t2=period-t1;
       DC=t1/period;
       par_grat.itetm=3;
       par_grat.px=period;
       par_grat.DC=DC;
       NModi=11;
       par_grat.NModi=NModi;
       itetm=3;
       icarico=0;   % calcola T reticolo caricato sui G_i
       dret=L_i(fish);
       %' cont grat', keyboard
       ibast=fish;
     else 
       ibast=[];
       dret=0;
     end  
   kt=0;
   

      if n_i(end)~=rfd
       L_i=[L_i; 0];
       n_i=[n_i; rfd];
      end   
   
   L_imic=L_i;
   

   
   dlam=.05;   % fraction of lam0 where res is searched
   [ga,la,zet,nz,Ez,Hz,]=CMM_1D(lam0,L_imic,n_i,iat,icav,ifp,dlam,rr,kt,ibast,par_grat,ipcam);