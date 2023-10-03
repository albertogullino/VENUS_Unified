function [Te,Tm,Nx,neq]=Tgrating(lambda,thick,rr,KK,par_grat,Nx_in)

segem=1;
iplotCa=0;

  if ~exist('Nx_in')
   Nx=1;
  else
   Nx=Nx_in;
  end
  neq=0;
  r_in=par_grat.r_in;
  r_out=par_grat.r_out;
  r1=par_grat.r1;
  r2=par_grat.r2;
  period=par_grat.per;
  DC=par_grat.DC;
  NModi=par_grat.NModi;
  d1=period*DC;
  d2=period*(1-DC);
  iret_BW=0;
  if isfield(par_grat,'iret_BW')==1
   iret_BW=par_grat.iret_BW;
  end
  Li=thick;
  
  if iret_BW==0
  Nx=1;
  phii=0;
  tetai=0;
  [T11,T22,T12,T21,s11]=orta_skewTOTr(phii,tetai,r_in,r_out,r1,r2,d1,d2,thick,lambda,NModi,0,iplotCa,segem);

%' cont gra', keyboard

Ttote=[T11(1,1) T12(1,1); T21(1,1) T22(1,1)];    % dalla formula sin*cos
Ttotm=[T11(2,2) T12(2,2); T21(2,2) T22(2,2)];


[Gei,Gmi,Tei,Tmi]=gaperdm(KK,0,lambda,[],[],[],[],0,r_in,rr,0,[],[],rr);

[Geu,Gmu,Teu,Tmu]=gaperdm(KK,0,lambda,[],[],[],[],0,rr,rr,0,[],[],r_out);

Gti=[Gei];
o=ones(size(Gti));
oM=diag(o);
Gtu=[Geu];
gM=-diag(Gtu);
Vu=[oM gM; gM oM]/Teu;

gM=-diag(Gti);
Vi=[oM gM; gM oM]/Tei;

Te=Vu*Ttote*Vi;


Gti=[Gmi];
o=ones(size(Gti));
oM=diag(o);
Gtu=[Gmu];
gM=-diag(Gtu);
Vu=[oM gM; gM oM]/Tmu;

gM=-diag(Gti);
Vi=[oM gM; gM oM]/Tmi;

Tm=Vu*Ttotm*Vi;

else

%  TE
 kt=KK;
 n=rr;
 k0=2*pi/lambda;
 be=j*k0*n; 
     bb=conj(sqrt(1-kt.^2));
     ZEv=(1./bb);
      if ite==0
       ZEv=(bb);
     end
 er1=r1^2;
 er2=r2^2;
 tt=period;
 t1=d1;
 t2=d2;
 
 ex=tt*er1*er2/(t2*er1+t1*er2);
 ey=(t1*er1+t2*er2)/tt;
 
 ni=sqrt(ey);
 if ite==1
  neq=ni;
 end 
 del=(ni^2-n^2)/(2*n^2)*ZEv;
 M=[-(bb+del) -del; del (bb+del)];
  if ite==0  
   delz=(1-(n/ni)^2)/2*ZEv*(kt/bb)^2;
   Mz=[-delz delz; -delz delz];
   M=M+Mz;
  end
 dx=Li/Nx;
 Mv=dx*be*M;
 Te=expm(Mv); 

 ni=sqrt(ex);
 if ite==2
  neq=ni;
 end  
 del=(ni^2-n^2)/(2*n^2)*ZEv;
 M=[-(bb+del) -del; del (bb+del)];
  if ite==0  
   delz=(1-(n/ni)^2)/2*ZEv*(kt/bb)^2;
   Mz=[-delz delz; -delz delz];
   M=M+Mz;
  end
 dx=Li/Nx;
 Mv=dx*be*M;
 Tm=expm(Mv);
end
