%
% function [gain,delf,lambda,Eqw,Eout,xro,fian,Pa,Plot,Ppol]=...
%  VELM(fil_str,par_in,Dlam_mod,nK_dis,alim_in,...
%       iLP,ipolar,Ev_or_Od,nmasce,numodiacc,...
%       xroI,fimaxi,lfi_inp,r_pil,ro_ini,zeta,N,T,...
%       iplan,idyn,i2D,iany,ianys,fileeps,ifp,iPERD_vol);
%
% USCITE
%
% gain:     vettore guadagni di soglia in 1/s
% delf:     vettore Delta_f/f di soglia
% lambda:   lung. d'onda (micron) di risonanza della struttura planare
% Eqw:      campi nel QW
% Eout:     campi nella sez. di uscita
% xro:      vettore coord. radiale
% fian:     vettore coord. angolare
% Pa, Plot, Ppol:       Parametri modali di uscita
%
% INGRESSI
%
% fil_str:      file struttura
% par_in:       parametri come da file di struttura
% Dlam_mod:     intervallo ricerca modi in nm (inizio, fine, # punti)
% nK_dis:       # modi k trasversale
% alim_in:      k trasversale massimo normalizzato (valori tipici 0.1:0.25)
% iLP:          soluzione scalare (iLP=1) o vettoriale (iLP=0)
% ipolar:       solo per soluzione vettoriale (-1,1,2 entrambe)
% Ev_or_Od:     modi pari o dispari; modo fondamentale per Ev_or_Od='Even'
% nmasce:       max. numero modi calcolati per ogni tipo
% numodiacc:    numero modi azimutali accoppiati
%
% xroI:     vettore coord. radiale. Se numero intero, da numero punti di xro
% fimaxi:   max valore angolo (di solito 2*pig, ma anch pig per partic. simmetrie)
% lfi_inp:  vettore coord. angolare. Se numero intero, da numero punti di fian
% r_pil:    raggio ossido o iniezione
% ro_ini:    vettore coord. radiale N e T
% zeta:     vettore coord. longitudinale N e T
%
% iplan:    1:guadagno piatto in QW  0: forma data da N o da file
% idyn:     1: calcola Pa per dinamica
% iany:     1: include effetto elettro ottico
% ianys:    1: include effetto elasto ottico
% fileeps:  nome file .mat con dati guadagno QW
% ifp:      flag di controllo e visualizzazione.
%           -4(no stop), -10(1 stop), 0..3(stop crescenti)

function [gain,delf,lambda,Eqw,Eout,xro,fian,Pa,Ppt,Ppo,ord]=...
VELM(fil_str,par_in,Dlam_mod,nK_dis,alim_in,...
 iLP1,ipolar,Ev_or_Od,nmasce,numodiacc,...
 xroI,fimaxi,lfi_inp,r_pil,ro_ini,zeta,N_in,T,...
 iplan,idyn,i2D,iany,ianys,fileeps,ifp,iPERD_vol,ifiint);

' invelm', keyboard

format compact
format short e

global ilo
global IdeOo IdeOon pMu0u pMu0u1 lKA nk1max nures pMc sim0 numodi pMei iLP ldap


global del_n_ag ianti_gui Nref
sN=size(N_in);
sN0=size(N_in);
igau=0;
yiN=0;
yiT=0;

%' qui sN', keyboard
if min(sN)==1
 N=N_in;
 sN=1;
else
 if sN(2)==2
  sN=1;
  N=N_in(:,1);
  Nref=N_in(:,2);
 else
  sN=ismat(N_in);
  N=N_in;
  igau=4;
  yiN=[1:100]; %solo per renderlo un vettore
 end
end
%' invelm', keyboard



if isstruct(ro_ini)==1
 ro_in=ro_ini.N;
 ro_inT=ro_ini.T;
else
 ro_in=ro_ini;
 ro_inT=ro_ini;
end


' VELM', keyboard
if length(ianti_gui)==0
 ianti_gui=0;
end
iraff=1;
ivfre0=0;
nvar=fil_str;
Dla=0;

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
 if ifp~=-4
  keyboard
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
 if iLP==1
  mmvet=[0 1];
 else
  mmvet=[1 0];
  
  if length(numodiacc)==1
   if numodiacc==0
      mmvet=[1 2];
   end   
  else
   if numodiacc(2)==0
    mmvet=[1 2];
   end
  end 
 end
else
 mmvet=str2num(Ev_or_Od);
end

%'vel', keyboard
%'vel', keyboard

global ilossk exp_los kl perdk
if length(ilossk)==0
ilossk=[0 0 0];
perdk=0;
%emme='emme_mix';
%emme='emme_sav';
emme='emme_ult';
%emme='em_mixp';
exp_los=1;
kl=.1;
else
%perdk=.1;
 if sum(ilossk)>0
  emme='em_mixp';
 else
%  emme='emme_mix';
  emme='emme_ult';
 end
%exp_los=25;
%kl=.35;
end

% N and T discretization

ndz=12;
ndz=20;


nst2=5;
nst1=10;
nst0=5;

fmol=2;
ndz=3*fmol;


nst2=5*fmol;
nst1=20*fmol;
nst0=7*fmol;

%ndz=20;
%nst0=7;
%nst1=15;
%nst2=12;

iold=0;
if iold==1
 if abs(min(zeta)-max(zeta))<15
  zdis=linspace(min(zeta),max(zeta),ndz+1);
 else
  dz=abs(min(zeta)-max(zeta))/6;
  [du,izm]=max(T(1,:));
  if izm>length(zeta)/2
   zmi=zeta(izm)-dz;
   zma=zeta(end);
   zdis=[zeta(1) linspace(zmi,zma,ndz)];
  else
   zma=zeta(izm)+dz;
   zmi=zeta(1);
   zdis=[linspace(zmi,zma,ndz) zeta(end)];
  end
 end
 T1=T(1,:);

if ifp>0  & length(zeta)>1
 figure, plot(zeta,T1,zdis,spline(zeta,T1,zdis),'r.'), pausak

% ' zeta ', keyboard
end

 zedis=zdis;

if length(ro_in)>1
 if length(r_pil)==0
  r_pil=max(ro_in)/3;
 end
 farpi=1.5;
 rdis0=linspace(0,r_pil/2,nst0);
% rdis1=linspace(r_pil/2,r_pil*farpi,nst+1);
% rdis2=linspace(r_pil*farpi,max(ro_in),nst1+1);
 rdis1=linspace(r_pil/2,r_pil*farpi,nst1+1);
 rdis2=linspace(r_pil*farpi,max(ro_in),nst2+1);
 rdisdu=[rdis0(1:end-1) rdis1 rdis2(2:end) rdis2(end)+diff(rdis2(1:2))];
 rdis=rdisdu;
%'rdis'
%keyboard
end


%' qui T 0', keyboard
if length(T)~=1

%' qui T', keyboard
 T0=(T(1,:));
 rodisT=rdis;
 for iz=1:length(zdis)-1
%  iz
% for iz=1:length(zdis)
  fi=find(zeta>=zdis(iz) & zeta<zdis(iz+1));
%  pausak
  Tz(iz)=mean(T0(fi));
  if length(fi)>1
   mea=mean(T(:,fi).');
  else
   mea=(T(:,fi).');
  end
  for ir=1:length(rdis)-1
   fi=find(ro_in>=rdis(ir) & ro_in<rdis(ir+1));
   Tr(ir)=mean(mea(fi));
  end
%  Tr(ir+1)=0;
  Tdis(:,iz)=Tr.';
 end
 fi=find(isnan(Tdis)==1);
 Tdis(fi)=0;
 rdis=rdis(1:end-1);


 zp=[zdis(1:ndz) zdis(2:ndz+1)-1e-4];
 tp=[Tz Tz];
 [zpu,iso]=sort(zp);
 tpu=tp(iso);

if ifp>0
 figure, plot(zpu,tpu,zeta,T0), pausak
 figure, plot(rdis,Tdis,'.-'), pausak
 if imag(Tdis(1))~=0
  figure, plot(rdis,imag(Tdis),'.-'),
  title(' Variazione perdite ')
  pausak
 end
% pausak
end

%idis_temp=1;
idis_temp=0;
'idistem', keyboard

if idis_temp==1

T0r=real(T(1,:));
dT=max(T0r)-min(T0r);
%keyboard
sdT=dT/(length(zdis)-1);

icd=1;
clear Tr Td zd Tdis1 zp tp
t0=T0r(1);
ki=1;
 for k=2:length(T0)
  if abs(T0r(k)-t0)>=sdT
   Td(icd)=mean(T0(ki:k));
   zd(icd)=zeta(k);
   t0=T0r(k);
  mea=mean(T(:,ki:k).');
  for ir=1:length(rdis)-1
   fi=find(ro_in>=rdis(ir) & ro_in<rdis(ir+1));
   Tr(ir)=mean(mea(fi));
  end
  Tr(ir+1)=mea(end);
  Tdis1(:,icd)=Tr.';
   ki=k;
   icd=icd+1;
  end
 end
 Td(end+1)=mean(T0r(ki:end));

 zp=[[0 zd(1:end)] [zd(1:end)-1e-4 zeta(end) ]];
 tp=[Td Td];
 [zpu1,iso]=sort(zp);
 tpu1=tp(iso);

if ifp>0
 figure, plot(zpu,tpu,zeta,T0r,zpu1,tpu1)
 title(' giallo: discr cost,  cian: variabile ')
 pausak
end


zedis=zd;
zdis=zd;

Tdis=[Tdis1; Tdis1(end,:)];

end  %idis_temp

rodisT=rdis(1:end-1);
% keyboard

 igau=4;
else
 rodisT=0;
 zdis=0;
 Tdis=0;
end

%' igau Temp ' , keyboard
global inuG

if inuG==1

igainshape=0;
yiN=0;
if sN==1

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
%  'passo N?', keyboard
   N=N/max(N);

   fir=find(rdis<=r_pil);
   %  for ir=1:length(fir)-1
    for ir=1:length(rdis)-1
     fi=find(ro_in>=rdis(ir) & ro_in<rdis(ir+1));

     rodisN(ir)=mean(rdis(ir:ir+1));
     if length(Nref)>1
      if ianti_gui==1
       Ndis_ag(ir)=mean(Nref(fi));
      end
     end
     if length(N)~=1
      Ndis(ir)=mean(N(fi));
     end
    end
    rodisN(ir+1)=rdis(ir+1);
    if length(N)~=1
     Ndis(ir+1)=0;
    end
     if length(Nref)>1
      if ianti_gui==1
       Ndis_ag(ir+1)=0;
      end
     end

         yiN=[];
         if length(Ndis)>1
          xiN=rodisN(2:length(rodisN));
          NN=Ndis/max(Ndis);
          yiN=-diff(NN);
         end
         if ianti_gui==1
          if length(Ndis_ag)>1
           xiN=rodisN(2:length(rodisN));
            NN_ag=Ndis_ag/max(abs(Ndis_ag));
%            NN_ag=Ndis_ag;
            yiN_ag=-diff(NN_ag);
          end
         end
if ifp==-10
 figure, plot(xiN,N(1)-cumsum(yiN),'.',ro_in,N,'-'), pausak
 if exist('yiN_ag')
%  figure, plot(xiN,Nref(1)-cumsum(yiN_ag),'.',ro_in,Nref,'-'), pausak
  Nrp=Nref/max(abs(Nref));
  figure, plot(xiN,Nrp(1)-cumsum(yiN_ag),'.',ro_in,Nrp,'-'), pausak
 end
end

 end

%' end N ', keyboard

 if T==0
  igau=-4;
 end

else
 Ndis=0;
 Ndis_ag=0;
 rodisN=0;
end

end %sN
% end N and T discretization
else  %nuovo N T

if length(r_pil)==0
 r_pil=max(ro_in)/3;
end

end %inuG
% N and T discretization

if ~exist('ndz')
 ndz=10;
end

famol=2;
nst=10*famol;
nst1=4*famol;
nst0=4*famol;

%nst=20;
%nst1=10;
%nst0=10;
 lT=ismat(T);
 if lT==3
  zdisd=linspace(min(zeta),max(zeta),ndz+2);
  zdis=(zdisd(1:end-1)+zdisd(2:end))/2;
 else
  zdis=linspace(min(zeta),max(zeta),ndz+2);
 end
 farpi=1.01;
 farpi=1.5;
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
  %keyboard
  pausak
  rdis=sort(rdu);
 end
 end


if sT>1

if ifp==-10
%' discretizzo Temp '
%keyboard

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
%  iz
%  pausak

  for ir=1:length(rdis)-1
   fi=find(ro_inT>=rdis(ir) & ro_inT<rdis(ir+1));
%   fi=find(ro_inT>=rdis(ir) & ro_inT<=rdis(ir+1));
%   fi=find(ro_inT>rdis(ir) & ro_inT<rdis(ir+1));
   Tr(ir)=mean(mea(fi));
  end

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
 figure, plot(rdis,Tdis,'.-')
% keyboard
pausak
 end
% ' qui dentro sT=2', keyboard
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
 zedis=(zdis(1:end-1)+zdis(2:end))/2;
end

'qui Temp', keyboard

igainshape=0;
sNt=[];
if sN==1
%' qui enne ', keyboard
yiN=0;
sN=ismat(N);
sNr=ismat(Nref);
%' sNr ', keyboard
sNt=[sN sNr];
if sN>0 | sNr>0
    if length(find(sNt)>=1)>=1
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
                    if length(Nref)~=1
                        if length(fi)>0
                            du1=mean(Nref(fi));
                        end
                        Ndis_ag(ir)=du1;
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
                    Ndis_ag(ir+1)=0;
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
%                    NN_ag=Ndis_ag/max(Ndis_ag);
                    yiN_ag=-diff(NN_ag);
                end
                if ifp==-10
                   if sN>0
                    figure, plot(xiN,N(1)-cumsum(yiN),'.',ro_in,N,'-'), pausak
                   end
                   if sNr>0
                    figure,
                    plot(xiN,max(Nref)*(Nref(1)/max(Nref)-cumsum(yiN_ag)),...
                           '.',ro_in,Nref,'-'), pausak
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


end  %sN
end  %sN

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
%emme='emme_mix';


disp(' inizio mvguu ')
%keyboard
%mvgu
sT=ismat(T);
mvguu
disp(' fine VELM ')
%keyboard
