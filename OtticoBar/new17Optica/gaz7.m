function [gazk,Iku,Ikb,fak,Tu,Tb,Gu,Gb,Lf,Lcav,ztot,Ez,Hz,indz,nmean,Perd]=...
         gaz7(fiQWi,fiCavi,L_i_nm,n_i,rr,nb,nu,lambda0,fr,kv,iLP,ifp,iff)
if length(fiCavi)==2
 fiCavi=fiCavi(1):fiCavi(end);
end
%ifp=2
%keyboard
%iff=1
%Lv=flipud(L_i);
%rv=flipud(n_i);
% prova plot
L_i=L_i_nm*1e-3;
if iff==1
 figure, plot(cumsum(L_i(fiCavi)),n_i(fiCavi),'.'), pausak
 keyboard
end
ru=nu;
rb=nb;
Lu=[];
Lb=[];
%ru=[];
%rb=[];
fiQW=fiQWi+length(Lu);
fiCav=fiCavi+length(Lu);
%Lv=[Lu; L_i; Lb];
Lv=L_i;
Lt=sum(Lv);
%rv=[ru; n_i; rb];
rv=n_i;
z=[0; (cumsum(Lv))];
NQW=length(fiQW);
iqw=round(NQW/2);
nstu=[1:fiQW(iqw)];
nstb=[1+fiQW(iqw):length(Lv)];
%
Lvb=flipud(Lv(nstb));
rvb=flipud(rv(nstb));
Lvu=(Lv(nstu));
rvu=(rv(nstu));
if length(fiQW)>1
 fiQWd=[1:2:length(fiQW)*2-1];
 fiQWu=length(Lvu)+1-fiQWd(1:iqw);
else
 fiQWd=[];
 fiQWu=length(Lvu);
end

rqw=rvu(fiQWu(1));

subfiqw=find(fiCav==fiQW(iqw));
ficavu=fiCav(1:subfiqw);

fiQWb=length(Lvb)-fiQWd(1:iqw-1);
subdu=length(fiCav)-subfiqw-1;
ficavb=(length(Lvb)-subdu):length(Lvb);

%if length(fiQWu)>length(fiQWd)
% fiQWb=length(Lvb)-fiQWd(1:iqw-1);
% subdu=length(fiCav)-subfiqw-1;
% ficavb=(length(Lvb)-subdu):length(Lvb);
%else
% ficavb=[];
% fiQWb=[];
%end
%' gaz7 ', keyboard

fiQW0=nstu(length(nstu));
%
%disp('gaz6 entro'), keyboard
%
zu=[0; (cumsum(Lvu))];
zb=[0; (cumsum(Lvb))];
mu=rr/nu;
mub=rr/nb;

j=sqrt(-1.);
flp=1-iLP;
k0=2*pi/lambda0;


sginu=1;

f1u=[1; sginu/mu];
f1b=[1; sginu/mub];
ic=1;
zi=0;
nz(1)=rv(1);
NPXi=19;

%' in gaz7', keyboard

f1o=[1; sginu/mu];
zio=0;
nzo=[];
f1d=[1; sginu/mub];
zid=0;
nzd=[];

for nk=1:length(kv)
 k=kv(nk);
 Itu=0;
 Iau=0;
 Icau=0;
 Icauh=0;
 ic=1;
 for nst=1:length(Lvu)
 %nst, pausak
  ic=ic+1;
  r=rvu(nst);
  mi=rr/r;
  sq=sqrt(1-mi^2.*k^2)';
  sqz=sqrt(1-flp*mi^2.*k^2)';
  fi=r/real(r);
  sqfi=sq*fi;
  dx=Lvu(nst);
%  pausak
  if iff==1
%function [fz,hz,z,nz]=fieldz(L_i,n_i,NPXi,rr,nb,nu,lambda0,fr,k,iLP)
   [f1o,zio,nzo]=fieldx(dx,r,NPXi,rr,sqfi,sqz,lambda0,fr,f1o,zio,nzo);
  end
%   nst, pausak
  be=k0*(1+fr)*r*sq*fi;
  dte=be*dx;
  KKi=[cos(dte),              j*sin(dte).*mi/sqz;...
       j*sin(dte)./mi*sqz,     cos(dte)];
%  if dx==100
%   [dx r]
%   keyboard
%  end

   nz(ic)=r;
%   f0(:,ic)=KKi*f0(:,ic-1);
   s2=sin(2*dte);
   Is=.5*[(dte+s2/2) (mi/sqz).^2*(dte-s2/2)];
%   Int=abs((Is*abs(f1u(:,ic-1)).^2/mi^2)/be);
%   Intzu(nst,1)=Int;
   Int=abs((Is*abs(f1u(:,ic-1)).^2/mi^2)/be);
   Intzu(nst,1)=Int;
   Intzup(nst)=Int*abs(mi^2)*imag(1/mi^2);
   Itu=Itu+Int;
   Is=.5*[ (sqz/mi).^2*(dte-s2/2) (dte+s2/2)];
   Inth=abs((Is*abs(f1u(:,ic-1)).^2)/be);
   Intzu(nst,2)=Inth;
   f1u(:,ic)=KKi*f1u(:,ic-1);
 end  %L_v
 icu=ic;
 miu=mi;

 ic=1;
 for nst=1:length(Lvb)
  ic=ic+1;
  r=rvb(nst);
  mi=rr/r;
  sq=sqrt(1-mi^2.*k^2)';
  sqz=sqrt(1-flp*mi^2.*k^2)';
  fi=r/real(r);
  dx=Lvb(nst);
  if iff==1
   [f1d,zid,nzd]=fieldx(dx,r,NPXi,rr,sqfi,sqz,lambda0,fr,f1d,zid,nzd);
  end
  be=k0*(1+fr)*r*sq*fi;
  dte=be*dx;
  KKi=[cos(dte),              j*sin(dte).*mi/sqz;...
       j*sin(dte)./mi*sqz,     cos(dte)];
   nz(ic)=r;
%   f0(:,ic)=KKi*f0(:,ic-1);
   s2=sin(2*dte);
   Is=.5*[(dte+s2/2) (mi/sqz).^2*(dte-s2/2)];
   Int=abs((Is*abs(f1b(:,ic-1)).^2/mi^2)/be);
   Intzb(nst,1)=Int;
   Intzbp(nst)=Int*abs(mi^2)*imag(1/mi^2);
   Is=.5*[ (sqz/mi).^2*(dte-s2/2) (dte+s2/2)];
   Inth=abs((Is*abs(f1b(:,ic-1)).^2)/be);
   Intzb(nst,2)=Inth;
   f1b(:,ic)=KKi*f1b(:,ic-1);
 end  %L_v
 icb=ic;
 mib=mi;
 fn=max(f1u(1,:))/max(f1b(1,:));
 fn=(f1u(1,icu))/(f1b(1,icb));
 fn2=abs(fn)^2;
 fnu=fn;
if iff==1
 lfo=length(f1o);
 lfd=length(f1d);
 fnc=f1o(1,lfo)/f1d(1,lfd);
 Ez=[f1o(1,:) fliplr(f1d(1,:))*fnc]/f1o(1,lfo);
 Hz=[-f1o(2,:) fliplr(f1d(2,:))*fnc]/f1o(1,lfo);
 indz=[nzo fliplr(nzd)];
 zd=zid;
 pu1=[2:lfd 1];
 zd=zd(pu1);
 ztot=cumsum([zio; flipud(zd)]);
figure, plot(ztot,abs(Ez).^2*3,ztot,abs(Hz).^2*3,ztot,indz)
%keyboard
else
 ztot=[];
 indz=[];
 Ez=[];
 Hz=[];
end

 Intzb=Intzb*fn2;
 Intzbp=Intzbp*fn2;


 Lca=sum(Lvu(ficavu))+sum(Lvb(ficavb));
 Icae=sum(Intzu(ficavu,1))+sum(Intzb(ficavb,1));
 Icah=sum(Intzu(ficavu,2))+sum(Intzb(ficavb,2));
if length(fiQWb)>0
 Iae=sum(Intzu(fiQWu,1))+sum(Intzb(fiQWb,1));
 Iah=sum(Intzu(fiQWu,2))+sum(Intzb(fiQWb,2));
else
 Iae=sum(Intzu(fiQWu,1));
 Iah=sum(Intzu(fiQWu,2));
% Lca=sum(Lvu(ficavu));
% Icae=sum(Intzu(ficavu,1));
% Icah=sum(Intzu(ficavu,2));
end


 Ica=(Icae+Icah)/2;
% Ica=Icae;

 Ite=sum(Intzu(:,1))+sum(Intzb(:,1));
 Itep(1)=sum(Intzup);
 Itep(2)=sum(Intzbp);
 Ith=sum(Intzu(:,2))+sum(Intzb(:,2));

nmean=real((sum(Intzu(:,1).*rvu)+sum(Intzb(:,1).*rvb))/Ite);

%nmean=real((sum(Intzu(ficavu,1).*rvu(ficavu))+sum(Intzb(ficavb,1).*rvb(ficavb)))/Ica);
%' nmean ', keyboard

 It=(Ite+Ith)/2;
% It=Ite;
 Perd=Itep/It;

 miqw2=(rr/rvu(fiQW0))^2;
 miqw2=1;
% keyboard

 Iau0=Intzu(fiQW0,1)*miqw2;
% Iau0=Intzu(fiQW0,1);
% Ia=(Iae+Iah)/2;
 Ia=Iae*miqw2;
 gazk0(nk,1)=Iau0/It;

 gazk(nk,1)=Ia/It;
 fak(nk,1)=Ia/(Iau0*NQW);

%' fak ', keyboard

 Ikb(nk,1)=It/abs(fnu)^2*abs(mub);
 Iku(nk,1)=It*abs(mu);

 miqw=rr/rqw;
 A0uf=(f1u(1,1)+f1u(2,1)*mu)/sqrt(mu)/2;      % b_u
 A0ub=(f1u(1,1)-f1u(2,1)*mu)/sqrt(mu)/2;      % a_u

 A0bf=(f1b(1,1)+f1b(2,1)*mub)/sqrt(mub)/2;
 A0bb=(f1b(1,1)-f1b(2,1)*mub)/sqrt(mub)/2;

 Auf=(f1u(1,icu-1)+f1u(2,icu-1)*miqw)/sqrt(miqw)/2; % b_u
 Aub=(f1u(1,icu-1)-f1u(2,icu-1)*miqw)/sqrt(miqw)/2; % a_u

 Abf=(f1b(1,icb)+f1b(2,icb)*miqw)/sqrt(miqw)/2;
 Abb=(f1b(1,icb)-f1b(2,icb)*miqw)/sqrt(miqw)/2;

 Tu(nk,1)=A0uf/Auf;
 Gu(nk,1)=Aub/Auf;
 Tb(nk,1)=A0bf/Abf;
 Gb(nk,1)=Abb/Abf;
 Lf(nk,1)=It/Ica*Lca;

Intud=[0 0; Intzu];
Intbd=[Intzb; 0 0];

% figure, plot(zu,real(f1u(1,:)),Lt-zb,real(f1b(1,:)*fn))
% figure, plot(zu,Intud,'.',Lt-zb,Intbd,'.')
% figure, plot(zu,f1u,'.',zb,f1b)
% disp(' fine k')
% pausak
end    %k
%pausak

% [[1:length(Lv)]' Lv rv]

%%fz(1,:)=(f0(1,:))/max(f0(1,:));
%fz(2,:)=(f1(1,:))/max(f1(1,:));
%%hz(1,:)=(f0(2,:))/max(f0(2,:));
%hz(2,:)=(f1(2,:))/max(f1(2,:));

%disp(' fieldz')
% figure, plot(cumsum(zi),abs(f0)), pausak
% figure, plot(z,abs(fz)), pausak
if ifp>0
 disp(' verifica gaz')
 vg=3e14/rr;
 peru1=abs(Tu)^2*vg/(2*Lf);
 peru2=vg/Iku;
 pero1=abs(Tb)^2*vg/(2*Lf);
 pero2=vg/Ikb;
 disp(' perdite uscita (con S21 e con definizione) ')
 [peru1 peru2]
 disp(' perdite sotto (con S21 e con definizione) ')
 [pero1 pero2]
 disp('gaz7')
 keyboard

end

%Lca=max(L_i);
%fi1=find(L_i==Lca);
%ficav=fi1(1):fi1(2);
Lcav=sum(L_i(fiCavi));

%disp('gaz7'), keyboard

if ifp>-3
% keyboard
end
