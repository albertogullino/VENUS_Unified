function [Eqw,Eout,xro,fian,lambda,delf,gain,ord,nrAzim,Cu,PaDy,...
          ADom,ADcm,Plot,Ppol,Dla_out,imod_err]=...
vel03(N,T,ro_in,zeta,xroI,fimaxi,lfi_inp,par_in,r_pil,nvar,...
 ipolar,Ev_or_Od,nmasce,numodiacc,Dlam_mod,nK_dis,alim_in,...
 iplan,iraff,idyn,i2D,iLP1,iany,ianys,ifp,Dla,fileeps,nomes,ivfre0,icampi);


%' vel', pausak
%   ipolar,Ev_or_Od,numodiacc,Dlam_mod,nK_dis,

% ipolar,Ev_or_Od,numodiacc,Frisi,Frisu,alim,Ndisp,nk1max1,...

format compact
format short e


global ilo ldapu
global IdeOo IdeOon pMu0u pMu0u1 lKA nk1max nures pMc sim0 numodi pMei iLP ldap

global del_n_ag ianti_gui Nref
if length(ianti_gui)==0
 ianti_gui=0;
end

 iLP=iLP1;
global Temper_sav 
Temper_sav=T;
%'set Temper_sav ', keyboard
%% Costanti universali

mass0=9.1e-31;
h=6.626e-34;
hbar=h/(2.*pi);
kB=1.38e-23;
j=sqrt(-1);
mi=pi*4e-7;
eps0=8.8541e-12;
c=1/sqrt(mi*eps0);
Z0=120*pi;
q=1.6e-19;
global pasnu

if length(pasnu)==0
 pasnu=2;
end
if pasnu==1
 ' !pasnu =1 '
 if ifp~=-4
%  keyboard
 end
end


if isequal(Ev_or_Od,'Even')
 if iLP==1
  mmvet=0;
 else
  mmvet=1;
 end
elseif isequal(Ev_or_Od,'Odd')
 if iLP==1
  mmvet=1;
 else
  mmvet=0;
 end
elseif isequal(Ev_or_Od,'Both')
 if iLP==0
  mmvet=[1 0];
 else
  mmvet=[0 1];
 end
else
 mmvet=0:str2num(Ev_or_Od);
end

%'vel'
%keyboard


% N and T discretization

ndz=1;

nst=10;
nst1=5;
nst0=3;

nst=20;
nst1=10;
nst0=10;

 zdis=linspace(min(zeta),max(zeta),ndz+1);
 zedis=(zdis(1:end-1)+zdis(2:end))/2;
 ro_inT=ro_in;
 farpi=1.01;
 farpi=2;
 rdis0=linspace(0,r_pil/2,nst0);
 rdis1=linspace(r_pil/2,r_pil*farpi,nst+1);
 rdis2=linspace(r_pil*farpi,max(ro_in),nst1+1);
 rdis=[rdis0(1:end-1) rdis1 rdis2(2:nst1+1)];
%'rdis'
%keyboard

igau=0;
if length(T)~=1

%' discretizzo Temp '
%keyboard
 T0=T(1,:);
 rodisT=rdis;
 for iz=1:length(zdis)-1
  fi=find(zeta>=zdis(iz) & zeta<zdis(iz+1));
  Tz(iz)=mean(T0(fi));
  mea=mean(T(:,fi)');
  for ir=1:length(rdis)-1
   fi=find(ro_in>=rdis(ir) & ro_in<rdis(ir+1));
   Tr(ir)=mean(mea(fi));
  end
  Tr(ir+1)=0;
  Tdis(:,iz)=Tr';
 end
 fi=find(isnan(Tdis)==1);
 Tdis(fi)=0;


 zp=[zdis(1:ndz) zdis(2:ndz+1)-1e-4];
 tp=[Tz Tz];
 [zpu,iso]=sort(zp);
 tpu=tp(iso);

 figure, plot(zpu,tpu,zeta,T0)
 figure, plot(rdis,Tdis,'.-')
% keyboard
 sT=ismat(T);

 igau=5;
else
 rodisT=0;
 zdis=0;
 Tdis=0;
end

igainshape=0;

if length(N)>1 | length(Nref)>1
 Ndis=0;
 Ndis_ag=0;

 if length(N)==2 & length(Nref)==2
  yiN=N;
  yiN_ag=Nref;
  xiN=ro_in;


%         yiN=[];
%         if length(Ndis)>1
%          xiN=rodisN(2:length(rodisN));
%          NN=Ndis/max(Ndis);
%          yiN=-diff(NN);
%         end
%         if ianti_gui==1
%          if length(Ndis_ag)>1
%           xiN=rodisN(2:length(rodisN));
%            NN_ag=Ndis_ag/max(Ndis_ag);
%            yiN_ag=-diff(NN_ag);
%          end
%         end


 else


   fir=find(rdis<=r_pil);
   %  for ir=1:length(fir)-1
    for ir=1:length(rdis)-1
     fi=find(ro_in>=rdis(ir) & ro_in<rdis(ir+1));
%     pausak

     rodisN(ir)=mean(rdis(ir:ir+1));
     if length(Nref)~=1
%      if ianti_gui==1
       if length(fi)>0
        du1=mean(Nref(fi));
       end
       Ndis_ag(ir)=du1;
%      end
     end
     if length(N)~=1
      if length(fi)>0
       du2=mean(N(fi));
      end
      Ndis(ir)=du2;
     end
    end
    rodisN(ir+1)=rdis(ir+1);
    if length(N)~=1
     Ndis(ir+1)=0;
    end
     if length(Nref)~=1
%      if ianti_gui==1
       Ndis_ag(ir+1)=0;
%      end
     end

         yiN=[];
         if length(Ndis)>1
          xiN=rodisN(2:length(rodisN));
          NN=Ndis/max(Ndis);
          yiN=-diff(NN);
         end
%         if ianti_gui==1
          if length(Ndis_ag)>1
           xiN=rodisN(2:length(rodisN));
            NN_ag=Ndis_ag/max(Ndis_ag);
            yiN_ag=-diff(NN_ag);
          end
%         end
if ifp==-10
 figure, plot(xiN,N(1)-cumsum(yiN),'.',ro_in,N,'-'), pausak
 figure, plot(xiN,Nref(1)-cumsum(yiN_ag),'.',ro_in,Nref,'-'), pausak
end

 end

%' end N '
%keyboard
 if T==0
  igau=-4;
 end

else
 Ndis=0;
 Ndis_ag=0;
 rodisN=0;
end

% end N and T discretization

nT=1;
Tvet=300;

Dla0=Dla;

itetm=0;  %=1 TE, =2 TM con mu=0;
if iLP==0
 itetm=0;
end

iriga=0;
ifnm=1;

global ivfre
if exist('ivfre0')
 ivfre=ivfre0;
else
 ivfre=1;
end

isi=0;
ikiaut=0;
emme='emme_mix';
emme='emme_ult';

if iLP==1
 nomes=[nomes,'LP'];
end

%disp(' inizio mvg_chaos '), keyboard
mvg_chaos
if ifp~=-4
disp(' fine mvg_chaos')
end

%if exist('Dla_new')
% Dlam_out=Dla_new;
%else
% Dlam_out=Dlam_mod;
%end