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
  if isfield(par_grat,'r1')==1
   r1=par_grat.r1;
   r2=par_grat.r2;
  else
   r1=par_grat.n1;
   r2=par_grat.n2;  
  end
  period=par_grat.per;
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
  if length(r1)>1
   Li=par_grat.th/1000;
  end
%' cont gra', keyboard  
  if iret_BW==0
  Nx=1;
  phii=0;
  tetai=0;
%kref=2*pi/lambda*rr;
%kout=2*pi/lambda*par_grat.r_out;
%kt=kout*sin(teta/180*pi)/kref;
tetai=asin(rr/par_grat.r_out*KK)*180/pi;

  [T11,T22,T12,T21,s11,s12,s21,s22]=orta_skewTOTr(phii,tetai,r_in,r_out,r1,r2,d1,d2,Li,lambda,NModi,0,iplotCa,segem);

%' cont gra', keyboard

Ttote=[T11(1,1) T12(1,1); T21(1,1) T22(1,1)];    % dalla formula sin*cos
Ttotm=[T11(2,2) T12(2,2); T21(2,2) T22(2,2)];
Stotm=[s11(2,2) s12(2,2); s21(2,2) s22(2,2)];

[Gei,Gmi,Tei,Tmi]=gaperdm(KK,0,lambda,[],[],[],[],0,rr,rr,0,[],[],r_in);

[Geu,Gmu,Teu,Tmu]=gaperdm(KK,0,lambda,[],[],[],[],0,r_out,rr,0,[],[],rr);

[Gei,Gmi,Tei,Tmi]=ga_simp(r_in,rr,KK,rr);

[Geu,Gmu,Teu,Tmu]=ga_simp(rr,r_out,KK,rr);

oM=eye(size(Gei));
gM=diag(Geu);
Vu=[oM gM; gM oM]/Teu;

gM=diag(Gei);
Vi=[oM gM; gM oM]/Tei;

Te=Vu*Ttote*Vi;
%Te=inv(Te1)';
%Te=Ttote;


gM=diag(Gmu);
Vu=[oM gM; gM oM]/Tmu;
gM=diag(Gmi);
Vi=[oM gM; gM oM]/Tmi;

Tm=Vu*Ttotm*Vi;
%Tm=inv(Tm1)';
%Tm=Ttotm;
%' qui Teq', keyboard

else

%  TE
 kt=KK;
 n=rr;
 k0=2*pi/lambda;
 be=j*k0*n; 
     bb=conj(sqrt(1-kt.^2));
     ZEv=(1./bb);
      if ite~=1
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
   if ite~=1  
   delz=(1-(n/ni)^2)/2*ZEv*(kt/bb)^2;
   Mz=[-delz delz; -delz delz];
   M=M+Mz;
  end
 dx=Li/Nx;
 Mv=dx*be*M;
 Te=expm(Mv); 

 ni=sqrt(ex);
 if ite==1
  neq=ni;
 end  
 del=(ni^2-n^2)/(2*n^2)*ZEv;
 M=[-(bb+del) -del; del (bb+del)];
  if ite~=1
   delz=(1-(n/ni)^2)/2*ZEv*(kt/bb)^2;
   Mz=[-delz delz; -delz delz];
   M=M+Mz;
  end
 dx=Li/Nx;
 Mv=dx*be*M;
 Tm=expm(Mv);
end


%' end Tgrating', keyboard