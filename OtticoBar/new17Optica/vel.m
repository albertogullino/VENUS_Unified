function [Eqw,Eout,xro,fian,lambda,delf,gain,ord,nrAzim,Cu,Pa,...
          ADom,ADcm,Plot,Ppol,Dla_out,imod_err]=...
vel_DD(Nref,N,T,ro_ini,zeta,xroI,fimaxi,lfi_inp,par_in,r_pil,nvar,...
 ipolar,Ev_or_Od,nmasce,numodiacc,Dlam_mod,nK_dis,alim_in,...
 iplan,iraff,idyn,i2D,iLP1,iany,ianys,ifp,Dla,fileeps,nomes,ivfre0,icampi,ndz,ivacc);

Eqw=0;
Eout=0;
xro=0;
fian=0;
lambda=0;
delf=0;
gain=0;
ord=0;
nrAzim=0;
Cu=0;
PaDy=0;
ADom=0;
ADcm=0;
Plot=0;
Ppol=0;
Dla_out=0;
imod_err=0;

global i1D

if ~exist('ivacc')
 ivacc=-1;
end
% iplan,iraff,idyn,i2D,iLP1,iany,ianys,ifp,Dla,fileeps,nomes,ivfre0,icampi,ro_inT);


%ifp=-10
%' vel', pausak
%end
if isstruct(ro_ini)==1
 ro_in=ro_ini.N;
 ro_inT=ro_ini.T;
else
 ro_in=ro_ini;
 ro_inT=ro_ini;
end
%' controllo', keyboard
%   ipolar,Ev_or_Od,numodiacc,Dlam_mod,nK_dis,

% ipolar,Ev_or_Od,numodiacc,Frisi,Frisu,alim,Ndisp,nk1max1,...

format compact
format short e

global ilo ldapu
global IdeOo IdeOon pMu0u pMu0u1 lKA nk1max nures pMc sim0 numodi pMei iLP ldap

global del_n_ag ianti_gui Nref ired_ret
if length(ianti_gui)==0
 ianti_gui=0;
end

 iLP=iLP1;

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
% if ifp~=-4
%  keyboard
% end
end
if iLP==0 & ipolar==0
% Ev_or_Od='Odd';
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
 mmvet=str2num(Ev_or_Od);
end

%'vel', keyboard


% N and T discretization

if ~exist('ndz')
 ndz=10;
end

famol=.5;
nst=10*famol;
nst1=6*famol;
nst0=8*famol;


'qui vel', keyboard
%nst=20;
%nst1=10;
%nst0=10;
 lT=ismat(T);
 if lT==3
  LzOpt=sum(L_i(1:end-1));
  zdisd=linspace(min(zeta),LzOpt,ndz+1);
  zdis=(zdisd(1:end-1)+zdisd(2:end))/2;
 else
  zdis=linspace(min(zeta),max(zeta),ndz+2);
 end
 farpi=1.01;
 farpi=max(ro_inT)/r_pil*.7;
 rdis0=linspace(0,r_pil/2,nst0);
 rdis1=linspace(r_pil/2,r_pil*farpi,nst+1);
 rdis2=linspace(r_pil*farpi,max(ro_inT),nst1+1);
 rdis=[rdis0(1:end-1) rdis1 rdis2(2:nst1+1)];
sT=ismat(T);
 if sT==2
 iprdis=0;
 if iprdis==1
  mt=mean(T,2);
  mt=mt(1:end-1)';
  xt=ro_inT(1:end);
  figure, plot(xt,mt)
  maT=max(mt);
  npr=length(rdis);
  dy=fix(100/npr)/100;
  dF=maT*dy/50;
  maT=maT-dF;
  miT=min(mt)+dF;
  mtp=linspace(miT,maT,npr);
  co=polyfit(mt,xt,30);
  rdu=polyval(co,mtp);
  figure, plot(xt,mt,rdu,mtp,'r.')

  'rdis'
  keyboard
  rdis=sort(rdu);
 end
 end

igau=0;

if sT>1

if ifp==-10
%' discretizzo Temp ', keyboard
end
 if sT==2

% ' qui dentro', keyboard
 T0=T(1,:);
 rodisT=rdis;
 for iz=1:length(zdis)-1
  fi=find(zeta>=zdis(iz) & zeta<zdis(iz+1));
%  fi=find(zeta>zdis(iz) & zeta<=zdis(iz+1));
  fiz=fi;
  Tz(iz)=mean(T0(fi));
  mea=mean(T(:,fi)');
%  iz, pausak

  for ir=1:length(rdis)-1
   fi=find(ro_inT>=rdis(ir) & ro_inT<rdis(ir+1));
%   fi=find(ro_inT>=rdis(ir) & ro_inT<=rdis(ir+1));
%   fi=find(ro_inT>rdis(ir) & ro_inT<rdis(ir+1));
   Tr(ir)=mean(mea(fi));
  end

%'qui Tr', keyboard
%  figure
% for iz=1:length(zdis)-1
%  fi=find(zeta>=zdis(iz) & zeta<zdis(iz+1));
%  fiz=fi;
%  Tz(iz)=mean(T0(fi));
%  mea=mean(T(:,fi)');
%  for ir=1:length(rdis)-1
%   fi=find(ro_inT>=rdis(ir) & ro_inT<rdis(ir+1));
%  end
% end

  Tr(ir+1)=mean(T(end,fiz));
  Tdis(:,iz)=Tr';
 end

%  figure
% for iz=1:length(zdis)-1
% [iz, zdis(iz)], pausak
%  fi=find(zeta>=zdis(iz) & zeta<zdis(iz+1));
%  fiz=fi;
%  Tz(iz)=mean(T0(fi));
%  mea=mean(T(:,fi)');
%  plot(ro_inT,T(1:end-1,fi),ro_inT,mea(1:end-1),'w.')
%  hold on
%  for ir=1:length(rdis)-1
%   fi=find(ro_inT>=rdis(ir) & ro_inT<rdis(ir+1));
%   Tri=mean(mea(fi));
%   plot(mean(ro_inT(fi)),Tri,'ro'), pausak
%  end
% end


 fi=find(isnan(Tdis)==1);
 Tdis(fi)=0;


 zp=[zdis(1:ndz+1) zdis(2:ndz+2)-1e-4];
 tp=[Tz Tz];
 [zpu,iso]=sort(zp);
 tpu=tp(iso);
 if ifp==-10
 figure, plot(zpu,tpu,zeta,T0), pausak
 figure, plot(rdis,Tdis,'.-'), pausak

% keyboard
 end
 %' qui dentro sT=2', keyboard
 elseif sT==3


 end

 igau=4;
else
 rodisT=0;
 zdis=0;
 Tdis=0;
end

if lT==3
 zedis=zdis;
else
% zedis=(zdis(1:end-1)+zdis(2:end))/2;
 zedisp=zdis(2:end);
 zedis=zdis;
end

%'qui zedis', keyboard

igainshape=0;

if ifp==-10
' qui enne ', keyboard
end
yiN=0;
sN=abs(ismat(N));
sNr=abs(ismat(Nref));
' N ', keyboard
if sN>0
    if min(sN)==1
        if length(N)>1 | length(Nref)>1
            Ndis=0;
            Ndis_ag=0;

            if length(N)==2 & length(Nref)==2
                yiN=N;
                yiN_ag=Nref;
                xiN=ro_in;

            else


                fir=find(rdis<=r_pil);
                %  for ir=1:length(fir)-1
                for ir=1:length(rdis)-1
                    fi=find(ro_in>=rdis(ir) & ro_in<rdis(ir+1));

                    rodisN(ir)=mean(rdis(ir:ir+1));
                    if length(Nref)>0
                    if length(Nref)~=1
                        if length(fi)>0
                            du1=mean(Nref(fi));
                        end
                        Ndis_ag(ir)=du1;
                    end
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
                if length(Nref)>0
                if length(Nref)~=1
                    Ndis_ag(ir+1)=0;
                end
                end
                yiN=[];
                if length(Ndis)>1
                    xiN=rodisN(2:length(rodisN));
                    NN=Ndis/max(Ndis);
                    yiN=-diff(NN);
                end
                if length(Ndis_ag)>1
                    xiN=rodisN(2:length(rodisN));
                    NN_ag=Ndis_ag/max(Ndis_ag);
                    yiN_ag=-diff(NN_ag);
                end
                if ifp==-10
                    figure, plot(xiN,N(1)-cumsum(yiN),'.',ro_in,N,'-'), pausak
                    if length(Nref)>0
                    figure, plot(xiN,Nref(1)-cumsum(yiN_ag),'.',ro_in,Nref,'-'), pausak
                    end
                end

%                ' end N '
%                keyboard
                if T==0
                    igau=-4;
                end
            end
        end
    else  % min sN
        yiN=[1:100]; %solo per renderlo un vettore
    end

else
    Ndis=0;
    Ndis_ag=0;
    rodisN=0;
end

%' dopo enne e T ', keyboard
 %end N and T discretization

nT=1;
Tvet=300;

Dla0=Dla;

itetm=-2;  %=1 TE, =2 TM con mu=0; neg per LP
%itetm=0;  %=1 TE, =2 TM con mu=0;
if iLP==0
 itetm=0;
end
%'itetm', keyboard
iriga=0;
ifnm=1;

global ivfre
if exist('ivfre0')
 ivfre=ivfre0;
else
 ivfre=1;
end

isi=0;
%isi=1;
% isi regola la simmetrizzazione del sistema, che si ottiene prendendo funzioni di base normalizzate
% a 1/sqrt(Cmu). Viene fatto dove si calcolano i coeff. di accoppiamento. Il programma 
% non e' ancora sistemato per i campi, che non tengono conto ancora della diversa base. Gli autovalori
% sono comunque uguali, cosa che verifica l'equivalenza.

% Altro parametro del sistema e' il dk. Questo viene regolato da iaut (in vel.m). Se iaut=0 si
% usa il metodo standard di includere il dk nella matrice di accoppiamento. Se iaut=1, invece 
% dk viene inglobato nell'autovettore (il che tuttavia comporta che dpes venga moltiplicato davanti
% alla matrice K. Anche in questo caso si provano i due metodi del tutto equivalenti.

%'set isi', keyboard
ikiaut=0;
global ilossk exp_los kl perdk
if length(ilossk)==0
ilossk=[0 0 0];
perdk=0;
%emme='emme_mix';
%emme='emme_sav';
%emme='emme_ult';
emme='emme_any';
emme='emme_anyB';
exp_los=1;
kl=.1;
else
%perdk=.1;
 if sum(ilossk)>0
  emme='em_mixp';
 else
%  emme='emme_mix';
%  emme='emme_ult';
 emme='emme_any';
 end
%exp_los=25;
%kl=.35;
end

if iLP==1
 nomes=[nomes,'LP'];
end

%emme='emme_sav';
%emme='emme_ult';
%emme='emme_any';
%emme='emme_anyB';
%emme='emme_anyB_old';
emme='emme_anyB_bot';
%emme='emme_any0';

global Pf

if isfield(Pf,'emme')==1
 emme=Pf.emme;
end

%disp(' inizio vel '), keyboard
%global Ps

mvg_last
if ifp~=-4
disp(' fine vel ')
end

%if exist('Dla_new')
% Dlam_out=Dla_new;
%else
% Dlam_out=Dlam_mod;
%end
