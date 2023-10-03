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

function [GGe,GGm,Te,Tm]=ga_grat(kv,fr,lambda0,Lv,rv,Lb,rb,n,rf,rr,...
                                Luv,ruv,ring,par_grat)

iLP=0;
dbstop if error

%1) m1=1.15537  m2=1.  ovvero no=3.495 n1=3.025 n2=n0   (in imped.m)

%3) k0L=pi  kiLi=pi/2  (negli elementi di matrice del metodo fogli)


j=sqrt(-1.);

mu=rr./rf;

if exist('ring')
 ming=rr/ring;
else
 ming=1;
end


flp=1-iLP;
%flp=1;
k0=2*pi/lambda0;
%' cont gae', keyboard
I=[1 0; 0 1];


for i=1:length(kv)
  k=kv(i);
 KKev=I;
 KKmv=I;

 for nst=1:length(Lv)
  Li=Lv(nst);
  r=rv(nst);
  mi=rr/r;
  sq=sqrt(1-mi'^2.*k^2)';
  sqz=sqrt(1-flp*mi'^2.*k^2)';
  pi1=k0*(1+fr)*Li*r;
  if Li>0
  KKie=[cos(pi1*sq),...
     j*sin(pi1*sq).*mi/sqz;...
     j*sin(pi1*sq)./mi*sqz,...
     cos(pi1*sq)];

   KKim=[cos(pi1.*sq),...
     j*sin(pi1.*sq).*mi*sqz;...
     j*sin(pi1.*sq)./mi/sqz,...
     cos(pi1.*sq)];
  else
  
  end
  KKev=KKev*KKie;
  KKmv=KKmv*KKim;

 end

 KKevu=I;
 KKmvu=I;
if exist('Luv') & length(Luv)>0
 for nst=1:length(Luv)
  Li=Luv(nst); 
  r=ruv(nst);
  mi=rr/r;
  sq=sqrt(1-mi'^2.*k^2)';
  sqz=sqrt(1-flp*mi'^2.*k^2)';
  pi1=k0*(1+fr)*Li*r;

  KKi=[cos(pi1*sq),...
     j*sin(pi1*sq).*mi/sqz;...
     j*sin(pi1*sq)./mi*sqz,...
     cos(pi1*sq)];
    if sq==0
     KKi(1,2)=j*mi/pi1;
    end     
%     [Luv(nst) r]
%     keyboard
  KKevu=KKevu*KKi;
%'qui e', keyboard
  if iLP==0
   KKi=[cos(pi1.*sq),...
     j*sin(pi1.*sq).*mi*sqz;...
     j*sin(pi1.*sq)./mi/sqz,...
     cos(pi1.*sq)];
    if sq==0
     KKi(2,1)=j/mi/pi1;
    end
%    'ver', keyboard
   KKmvu=KKmvu*KKi;
  end
 end
% 'qui m', keyboard
end

 KKe=I;
 KKm=I;
 if n~=0
   for nst=1:length(Lb)
    Li=Lb(nst);   
    r=rb(nst);
    mi=rr/r;
    sq=sqrt(1-mi'^2.*k^2)';
    sqz=sqrt(1-flp*mi'^2.*k^2)';
    pi1=k0*(1+fr)*Li*r;

    KKi=[cos(pi1*sq),...
       j*sin(pi1*sq).*mi/sqz;...
       j*sin(pi1*sq)./mi*sqz,...
       cos(pi1*sq)];
       
       
    if sq==0
     KKi(1,2)=j*mi/pi1;
    end

    KKe=KKe*KKi;



   if iLP==0
    KKi=[cos(pi1.*sq),...
       j*sin(pi1.*sq).*mi*sqz;...
       j*sin(pi1.*sq)./mi/sqz,...
       cos(pi1.*sq)];
    if sq==0
     KKi(2,1)=j/mi/pi1;
    end
    KKm=KKm*KKi;
   end
   
    if nst==1
     KKe1=KKe;
     KKm1=KKm;
    end
   end  %nst
 end  %n

     zue=(mu/sqrt(1-flp*k.^2.*mu'.^2))';
     z0e=(ming/sqrt(1-flp*k.^2*ming'^2))';
     zum=(mu*sqrt(1-flp*k.^2.*mu'.^2))';
     z0m=(ming*sqrt(1-flp*k.^2*ming^2))';

%     zue=(mu/sqrt(1-flp*k.^2.*mu'.^2));
%     z0e=(ming/sqrt(1-flp*k.^2*ming'^2));
%     zum=(mu*sqrt(1-flp*k.^2.*mu'.^2));
%     z0m=(ming*sqrt(1-flp*k.^2*ming^2));

%     beu=conj(sqrt(1-flp*k.^2.*mu^2));
%     bei=conj(sqrt(1-flp*k.^2*ming^2));
%     zue=mu/beu;
%     z0e=ming/bei;
%     zum=mu*beu;
%     z0m=ming*bei;

     nint=fix(n);
% prova potenza matrice
% pM=(KKe)^nint;
% pMf=[1 0; 0 1];
% for k=1:nint
%  pMf=(KKe)*pMf;
% end


     M1=KKev*((KKe)^nint);
     if n-nint~=0
      M1=M1*KKe1;
     end
     M1=M1*KKevu;
     A=(M1(1,1)*zue+M1(1,2))/z0e;
     B=(M1(2,1)*zue+M1(2,2));
     R=(A/B-1)/(A/B+1);
     GGe(i,1)=R;
     T=2*sqrt(real(zue/z0e))/(A+B);
     Te(i,1)=T;
%'quie', keyboard

%     A=(M1(1,1)*zue+M1(1,2))/z0e;
%     B=(M1(2,1)*zue+M1(2,2));
%     R=(A/B-1)/(A/B+1);
%     GGe(i,1)=R;
%     T=2*sqrt(real(z0e*zue))/(A*z0e+B*z0e);
%     Te(i,1)=T;
%  [abs(T)^2 abs(R)^2 1-abs(R)^2-abs(T)^2]
%  f0=[1; -1/zue];
%  fu=M1*f0
%     A=(M1(1,1)*zue-M1(1,2))/z0e;
%     B=(M1(2,1)*zue-M1(2,2));
%     R1=(A/B+1)/(A/B-1)
%     disp ('gaemms')
%     keyboard

%     M1=(KKe)^nint;
%     if n-nint~=0
%      M1=M1*KKe1;
%     end
%     A=(M1(1,1)*zue+M1(1,2))/z0e;
%     B=(M1(2,1)*zue+M1(2,2));
%     R=(A/B-1)/(A/B+1);

  if iLP==0
     M1=KKmv*(KKm)^nint;
     if n-nint~=0
      M1=M1*KKm1;
     end
     M1=M1*KKmvu;
     A=(M1(1,1)*zum+M1(1,2));
     B=(M1(2,1)*zum+M1(2,2));
     R=(A/B-z0m)/(A/B+z0m);
     GGm(i,1)=R;
     T=2*zum/(A+B*z0m);
%     T=2*sqrt(real(z0m*zum))/(A+B*z0m);
     Tm(i,1)=T;
   else
    Tm=Te;
    GGm=GGe;
   end
% k, pausak

 end
%if kv>.5
% 'gaemms new02',  keyboard
%  pausak
%end
