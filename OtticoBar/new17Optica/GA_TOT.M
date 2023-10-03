% function [GGe,GGm,Te,Tm]=
%   gaemms(kv,fr,lambda0,Lv,rv,Lb,rb,n,rf,rr,iLP,Luv,ruv,ring)
%
%  Lv, rv  : strati prima dello speccio
%  Lb, rb  : strati dello speccio
%  rf      : n_rif uscita
%  rr      : n_rif riferimento
%  Luv, ruv: strati dopo lo speccio
%
% gli strati vanno inseriti con la stessa  direzione verso cui si vuole
% calcolare Gamma.
% Es. Gamma visto da A degli strati A, b, c, D ---> ordine A, b, c, D
% Es. Gamma visto da D degli strati A, b, c, D ---> ordine D, c, b, A

function [GGe,GGm,Te,Tm,du,nu]=ga_tot(kv,fr,lambda0,Lvi,rvi,Lbi,rbi,ni,...
                                Luvi,ruvi,fape,rfi,rr,iLP)



dbstop if error

%1) m1=1.15537  m2=1.  ovvero no=3.495 n1=3.025 n2=n0   (in imped.m)

%3) k0L=pi  kiLi=pi/2  (negli elementi di matrice del metodo fogli)


j=sqrt(-1.);
%fape=(1-1./(1+(kv/kl).^18))*perd*(-j);

du=[];
nu=[];


flp=1-iLP;
%flp=1;
k0=2*pi/lambda0;

I=[1 0; 0 1];





for i=1:length(kv)

%sopra
Lv=Lvi.t;
rv=rvi.t;
Lb=Lbi.t;
rb=rbi.t;
Luv=Luvi.t;
ruv=ruvi.t;

n=ni.t;
ring=rfi.t;
 k=kv(i);
 KKev=I;
 KKmv=I;
if i==1
 du=[.2];
 nu=[ring];
end


 if exist('ring')
  ringp=ring+fape(i);
  ming=rr/ringp;
 else
  ming=1;
 end

 for nst=1:length(Lv)
  r=rv(nst)+fape(i);
  if i==1
  du=[du Lv(nst)];
  nu=[nu r];
  end

  mi=rr/r;
  sq=sqrt(1-mi'^2.*k^2)';
  sqz=sqrt(1-flp*mi'^2.*k^2)';
  pi1=k0*(1+fr)*Lv(nst)*r;

  KKi=[cos(pi1*sq),...
     j*sin(pi1*sq).*mi/sqz;...
     j*sin(pi1*sq)./mi*sqz,...
     cos(pi1*sq)];
  KKev=KKev*KKi;
  if iLP==0
   KKi=[cos(pi1.*sq),...
     j*sin(pi1.*sq).*mi*sqz;...
     j*sin(pi1.*sq)./mi/sqz,...
     cos(pi1.*sq)];
   KKmv=KKmv*KKi;
  end
 end

 KKevu=I;
 KKmvu=I;
if exist('Luv') & length(Luv)>0
 for nst=1:length(Luv)
%  r=ruv(nst);
  r=ruv(nst)+fape(i);
  if i==1
  du=[du Luv(nst)];
  nu=[nu r];
  end
  mi=rr/r;
  sq=sqrt(1-mi'^2.*k^2)';
  sqz=sqrt(1-flp*mi'^2.*k^2)';
  pi1=k0*(1+fr)*Luv(nst)*r;

  KKi=[cos(pi1*sq),...
     j*sin(pi1*sq).*mi/sqz;...
     j*sin(pi1*sq)./mi*sqz,...
     cos(pi1*sq)];
%     [Luv(nst) r]
%     keyboard
  KKevu=KKevu*KKi;

  if iLP==0
   KKi=[cos(pi1.*sq),...
     j*sin(pi1.*sq).*mi*sqz;...
     j*sin(pi1.*sq)./mi/sqz,...
     cos(pi1.*sq)];
   KKmvu=KKmvu*KKi;
  end
 end
end

 KKe=I;
 KKm=I;
 if n~=0
   for nst=1:length(Lb)
%    r=rb(nst);
    r=rb(nst)+fape(i);
  if i==1
    du=[du Lb(nst)];
    nu=[nu r];
  end
    mi=rr/r;
    sq=sqrt(1-mi'^2.*k^2)';
    sqz=sqrt(1-flp*mi'^2.*k^2)';
    pi1=k0*(1+fr)*Lb(nst)*r;

    KKi=[cos(pi1*sq),...
       j*sin(pi1*sq).*mi/sqz;...
       j*sin(pi1*sq)./mi*sqz,...
       cos(pi1*sq)];
    KKe=KKe*KKi;


   if iLP==0
    KKi=[cos(pi1.*sq),...
       j*sin(pi1.*sq).*mi*sqz;...
       j*sin(pi1.*sq)./mi/sqz,...
       cos(pi1.*sq)];
    KKm=KKm*KKi;
   end
    if nst==1
     KKe1=KKe;
     KKm1=KKm;
    end
   end  %nst
 end  %n


     nint=fix(n);
% prova potenza matrice
% pM=(KKe)^nint;
% pMf=[1 0; 0 1];
% for k=1:nint
%  pMf=(KKe)*pMf;
% end


      M1e=KKev*((KKe)^nint);
      if n-nint~=0
       M1e=M1e*KKe1;
      end
      M1e=M1e*KKevu;

     if iLP==0
      M1m=KKmv*((KKm)^nint);
      if n-nint~=0
       M1m=M1m*KKm1;
      end
      M1m=M1m*KKmvu;
     end


%sotto
Lv=Lvi.b;
rv=rvi.b;
Lb=Lbi.b;
rb=rbi.b;
Luv=Luvi.b;
ruv=ruvi.b;

n=ni.b;
rf=rfi.b;


 rfp=rf+fape(i);
 mu=rr./rfp;
 k=kv(i);
 KKev=I;
 KKmv=I;

 rfp=rf+fape(i);
 mu=rr./rfp;

 if exist('ring')
  ringp=ring+fape(i);
  ming=rr/ringp;
 else
  ming=1;
 end


 for nst=1:length(Lv)
  r=rv(nst)+fape(i);
  if i==1
  du=[du Lv(nst)];
  nu=[nu r];
  end
  mi=rr/r;
  sq=sqrt(1-mi'^2.*k^2)';
  sqz=sqrt(1-flp*mi'^2.*k^2)';
  pi1=k0*(1+fr)*Lv(nst)*r;

  KKi=[cos(pi1*sq),...
     j*sin(pi1*sq).*mi/sqz;...
     j*sin(pi1*sq)./mi*sqz,...
     cos(pi1*sq)];
  KKev=KKev*KKi;
  if iLP==0
   KKi=[cos(pi1.*sq),...
     j*sin(pi1.*sq).*mi*sqz;...
     j*sin(pi1.*sq)./mi/sqz,...
     cos(pi1.*sq)];
   KKmv=KKmv*KKi;
  end
 end

 KKe=I;
 KKm=I;
 if n~=0
   for nst=1:length(Lb)
%    r=rb(nst);
    r=rb(nst)+fape(i);
  if i==1
    du=[du Lb(nst)];
    nu=[nu r];
  end
    mi=rr/r;
    sq=sqrt(1-mi'^2.*k^2)';
    sqz=sqrt(1-flp*mi'^2.*k^2)';
    pi1=k0*(1+fr)*Lb(nst)*r;

    KKi=[cos(pi1*sq),...
       j*sin(pi1*sq).*mi/sqz;...
       j*sin(pi1*sq)./mi*sqz,...
       cos(pi1*sq)];
    KKe=KKe*KKi;


   if iLP==0
    KKi=[cos(pi1.*sq),...
       j*sin(pi1.*sq).*mi*sqz;...
       j*sin(pi1.*sq)./mi/sqz,...
       cos(pi1.*sq)];
    KKm=KKm*KKi;
   end
    if nst==1
     KKe1=KKe;
     KKm1=KKm;
    end
   end  %nst

 end  %n

 KKevu=I;
 KKmvu=I;
if exist('Luv') & length(Luv)>0
 for nst=1:length(Luv)
%  r=ruv(nst);
  r=ruv(nst)+fape(i);
  if i==1
    du=[du Luv(nst)];
    nu=[nu r];
  end
  mi=rr/r;
  sq=sqrt(1-mi'^2.*k^2)';
  sqz=sqrt(1-flp*mi'^2.*k^2)';
  pi1=k0*(1+fr)*Luv(nst)*r;

  KKi=[cos(pi1*sq),...
     j*sin(pi1*sq).*mi/sqz;...
     j*sin(pi1*sq)./mi*sqz,...
     cos(pi1*sq)];
%     [Luv(nst) r]
%     keyboard
  KKevu=KKevu*KKi;

  if iLP==0
   KKi=[cos(pi1.*sq),...
     j*sin(pi1.*sq).*mi*sqz;...
     j*sin(pi1.*sq)./mi/sqz,...
     cos(pi1.*sq)];
   KKmvu=KKmvu*KKi;
  end
 end
end


     zue=(mu/sqrt(1-flp*k.^2.*mu'.^2))';
     z0e=(ming/sqrt(1-flp*k.^2*ming'^2))';
     zum=(mu*sqrt(1-flp*k.^2.*mu'.^2))';
     z0m=(ming*sqrt(1-flp*k.^2*ming^2))';

     zue=(mu/sqrt(1-flp*k.^2.*mu'.^2));
     z0e=(ming/sqrt(1-flp*k.^2*ming'^2));
     zum=(mu*sqrt(1-flp*k.^2.*mu'.^2));
     z0m=(ming*sqrt(1-flp*k.^2*ming^2));

     nint=fix(n);
% prova potenza matrice
% pM=(KKe)^nint;
% pMf=[1 0; 0 1];
% for k=1:nint
%  pMf=(KKe)*pMf;
% end


     M1e=KKev*((KKe)^nint)*KKevu*M1e;

     A=(M1e(1,1)*zue+M1e(1,2))/z0e;
     B=(M1e(2,1)*zue+M1e(2,2));
     R=(A/B-1)/(A/B+1);
     GGe(i,1)=R;
     T=2*sqrt(real(z0e*zue))/(A*z0e+B*z0e);
     Te(i,1)=T;


  if iLP==0
     M1=KKmv*((KKm)^nint)*KKmvu*M1m;

     A=(M1(1,1)*zum+M1(1,2))/z0m;
     B=(M1(2,1)*zum+M1(2,2));
     R=(A/B-1)/(A/B+1);
     GGm(i,1)=R;
     T=2*sqrt(real(z0m*zum))/(A*z0m+B*z0m);
     Tm(i,1)=T;
   else
    Tm=Te;
    GGm=GGe;
   end
% k, pausak

 end
%if kv>.5
 du=[du .2];
 nu=[nu rf];
 cdu=cumsum(du);
 du=[[0 cdu(1:end-1)]; cdu];
 du=reshape(du,prod(size(du)),1);
 nu=[nu; nu];
 nu=reshape(nu,prod(size(nu)),1);
% figure, plot(du,nu,'r',du,nu,'w.')
% 'gaemms new02'
% keyboard
%  pausak
%end
