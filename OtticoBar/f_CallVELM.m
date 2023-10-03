
function [velm] = f_CallVELM(mesh,mode,mode1,ParVet,VelmOptions,fil_str)

clear Ps
clear global Ps

load saclearglo
eval(stringaclear)
%clear global

warning off
%==========================================================================
% VELM solver parameters
%==========================================================================
verbVELM=abs(mode.verbVELM);

fclose('all') 
num_azim=VelmOptions.num_azim;
ivett=VelmOptions.ivett;
ipolar=VelmOptions.ipolar;
itutmir=VelmOptions.itutmir;
if isfield(VelmOptions,'imod_acc')==1
 imod_acc=VelmOptions.imod_acc; 
else
 imod_acc=0; 
end

if isfield(VelmOptions,'ianti_gui')==1
 ianti_gui=VelmOptions.ianti_gui;
else 
 ianti_gui=1;
end

if isfield(VelmOptions,'Pf')==1
 Ps.Pf=VelmOptions.Pf;
end

if isfield(VelmOptions,'PaTJ')==1
 Ps.Pf.PaTJ=VelmOptions.PaTJ;
end


if isfield(VelmOptions,'dissfun')==1
 Ps.dissfun=VelmOptions.dissfun;
end
if isfield(VelmOptions,'igraef_new')==1
 Ps.igraef_new=VelmOptions.igraef_new;
end
if isfield(VelmOptions,'iBWnew')==1
 Ps.iBWnew=VelmOptions.iBWnew;
end
if isfield(VelmOptions,'iany')==1
 Ps.iany=VelmOptions.iany;
end

Ps.asspar=0;  %come prima
if mode.flgBTJ_lithographic>=1
 Ps.asspar=2;  %TJ vera
end


%'inati', keyboard
Ps.NP_sopra=15;
Ps.NP_sotto=2;

NPx=1e4*mesh.xgrid;

NPxQW=mesh.nnxQW{1};

dox=ParVet(1);

'NPx'%, keyboard
ired_int=0   % 1 riduce la zona dei fenomeni trasversali
% Fat_ox=2.5;
% Fat_ox=4;
Fat_ox=3;
if ired_int==1
%  Rmax=max([Fat_ox*dox 8]);
 Rmax=Fat_ox*dox;
 fianti=find(NPx<=Rmax);
else
 fianti=1:NPxQW;
 
 %% Here the full DD radial mesh is used instead of a portion: be careful to N and Nref!
%  fianti=1:length(mesh.xgrid);
end

% x=linspace(min(xdd),max(xdd),101);
% calop.ro_in=x';

if length(mesh.DeltaTvelm)>1
TddFl=reshape(mesh.DeltaTvelm,mesh.nny,mesh.nnx);
TddFl=TddFl(:,fianti);
else
TddFl=mesh.DeltaTvelm;
end
 Tdd=flipud(TddFl);
 
 T300=mode.T300; 

 DT_Peltier=mode.T0-T300;
 
% ATTENZIONE: tenere conto di dove si aumenta la T
Tdd=Tdd(:,fianti)+DT_Peltier;  % riduco la zona termica trasversale a dove il campo la sente
% Tdd=Tdd(:,fianti);  % riduco la zona termica trasversale a dove il campo la sente
  
 T0=max(max(Tdd));

%  Dlam(1)=Dlam(1)+T0/80;
%  Dlam(2)=Dlam(2)+T0/20;

% 'stop T0', keyboard
if T0<mode.DT0/4
% 'stop T0', keyboard
 %   if ParVet(1)<1
  %   'Ossido troppo piccolo: lo allargo a T=0 !!!!'
   %  ParVet(1)=1;
  % end
end


ro_in=NPx(fianti);

ro_Fie=NPx;

nnxVELM=351;    % original VELM for 2019JSTQE articles
% nnxVELM=401;
ro_campo=linspace(0,NPx(fianti(end)),nnxVELM);
ro_campoF=linspace(0,NPx(end),nnxVELM);



if length(mode.efield_y)>1
Es=reshape(mode.efield_y,mesh.nny,mesh.nnx);
else
Es=mode.efield_y;
end
Es=Es(:,fianti);



lambda=mean(mode.vlambda)*1e-9;

if length(mode.Deltan)>1
 DeN=mode.Deltan;
 if isfield(mode,'TARde')
  faDe=mode.TARde;
  ICAV=mesh.ICAV;
  DeN(ICAV)=faDe*DeN(ICAV);
  'ICAV'
 %keyboard
 else
  faDe=1;
 end

%    Deltan=reshape(mode.Deltan,mesh.nny,mesh.nnx);
    Deltan=reshape(DeN,mesh.nny,mesh.nnx);
    alphaDu=reshape(mode.alpha,mesh.nny,mesh.nnx);
    alphaDu=alphaDu(:,fianti);
    alpha_ax=alphaDu(:,1);
    fia=find(alpha_ax>1);
    alphaDue=zeros(size(alpha_ax));
    alphaDue(fia)=alphaDu(fia,end);
    alphaDueP=alphaDu(:,end);
    alphaEndP=repmat(alphaDueP,1,length(fianti));
    alphaEnd=repmat(alphaDue,1,length(fianti));
    alpha=alphaDu-alphaEnd;
    alphaP=alphaDu-alphaEndP;
    k0=2*pi/lambda;
    Deltan=Deltan(:,fianti);
%    ImagDeltan=alpha./(2.*k0).*(1+Tdd/T300).^mode.ABS_Texp;
    ImagDeltan=alpha./(2.*k0);
%'Deltan', keyboard
else
    Deltan=0;
    ImagDeltan=0;
end
%'Deltan', keyboard


Ex=reshape(mode.efield_x,mesh.nny,mesh.nnx); % rho field component
Ey=reshape(mode.efield_y,mesh.nny,mesh.nnx); % z field component
Ex=Ex(:,fianti);
Ey=Ey(:,fianti);

% for kz=1:length(calop.zeta)
%  Tvelm(:,kz)=spline(xdd,Tdd(kz,:),x);
% end
%calop.N=epsTf(:,1)/epsTf(1,1);

global inuT inuG
        inuT=1;
        inuG=1; 

NPx1=length(NPx);
iLP=1-ivett;
% fil_str='jim_ox1.str';

for kk=1:length(ParVet)
 par_in{kk,1}=ParVet(kk);
end
% parametro con tag |P1 nel file str

%'callvelm', keyboard
verbVELM=abs(mode.verbVELM);

if(verbVELM==2)
    ifp_Opt=-10;  % ifp_Opt=-10: modalità verbose, da usare con strutture nuove
else
    ifp_Opt=-4;   % ifp_Opt=-4: modalità veloce
end

if(verbVELM==0)
     Ps.ifpstop=0; 
end

%Dlam=[-0.2 1.3 7 0.3];   %scansione in frequenza per trovare il modo
%alim=.15;               %kmax; 0.1 va bene per aperture normali (3-4 um)
%NP_k=20;
Dlam=VelmOptions.Dlam;   %scansione in frequenza per trovare il modo
alim=VelmOptions.krel_max;               %kmax; 0.1 va bene per aperture normali (3-4 um)
NP_k=VelmOptions.NP_k;

% mod_max=[2];   % numero modi 
mod_max=[mode.nmodes];   % numero modi - Tibaldi
numodiacc=0;        %numero armoniche azimutali accoppiate
if isfield(VelmOptions,'numodiacc')==1
 numodiacc=VelmOptions.numodiacc;        %numero armoniche azimutali accoppiate
end
pasnu=2;             % passo armoniche, sempre 2 per strutture circolari

dk1=alim/NP_k;
% 'qui modi', keyboard
if ivett==1

if num_azim==2
 Ev_or_Od='Both';   % 'Even' se solo modo fondamentale, 'Odd' se solo modo superiore
elseif num_azim==1
 Ev_or_Od='Even';   % 'Even' se solo modo fondamentale, 'Odd' se solo modo superiore
else
 Ev_or_Od=num2str(num_azim-2); 
end

else
 if imod_acc==1
  pasnu=1;
  numodiacc=num_azim-1;
 end
 
 if num_azim==2
  Ev_or_Od='Both';   % 'Even' se solo modo fondamentale, 'Odd' se solo modo superiore
 elseif num_azim==1
  Ev_or_Od='Even';   % 'Even' se solo modo fondamentale, 'Odd' se solo modo superiore
 else
  Ev_or_Od=num2str(num_azim-1); 
 end 
end

% 'qui acc', keyboard
% sistemare
lambda=lambda*1e6;  % in um
rrat=3.6;

%%



%'g0', keyboard
%matDn=mode.matDn; % Tibaldi
%g0=abs(matgain(end));  % in 1/cm % Tibaldi
%g0=0;

setVELMparameters %-- chiamare i parametri che non devo settare

if VelmOptions.ianti_gui==1

  if isfield(mode,'DeltaN')==0
   Dn=0;
  else 
   Dn=mode.DeltaN;
   if length(Dn)>length(fianti)
    Dn=Dn(fianti);
   end 
  end
  
  del_n_agInf=Dn(end);
  del_n_ag0=Dn(1);
  del_n_ag=del_n_ag0-del_n_agInf;
  if del_n_ag==0
   Nref=0;
  else 
   Nref=(Dn-del_n_agInf)/del_n_ag;
   
   if fianti(end)>NPxQW
       'Nref=0'%,keyboard
       Nref(end+1:fianti(end))=0;
   end
   del_n_ag=del_n_ag*mode.fat_ag; 

   %del_n_ag=del_n_ag*(1-Nref_du(end));
%   'cont dfat[', keyboard
   
  end 

else
 Nref=0;
 del_n_ag=0;
 del_n_agInf=0;
end  


matgain=mode.matgain; % Tibaldi

if isfield(VelmOptions,'gain_gui')
 igain=VelmOptions.gain_gui;
 if igain==0
  matgain=0;
 end
end
if length(matgain)==1
    N=0;
%    vph=c_light*100/rrat;
%    g0=abs(matgain)/vph;
    g0=abs(matgain);
    iga=0;
    if iga==1
     roadd=8;
     ro_in=sort([ro_in roadd+.01]);
     fiN=find(ro_in<=roadd);
     N=zeros(size(ro_in));
     N(fiN)=1;
    else
     g0=0;
    end 
    %g0=0;
%    'quio', keyboard
else
    if length(matgain)>length(fianti)
        matgain=matgain(:,fianti);
    end
    if max(matgain)>-100
%         g0=abs(min(matgain));  % in 1/cm % Tibaldi
        N=matgain-min(matgain); N=N./max(N);
        if N(1)<0
            N=-N;
        end
        g0=abs((matgain(end)));  % in 1/cm % Tibaldi
        N=matgain-matgain(end); N=N./max(N);
        
        if fianti(end)>NPxQW
            'N=0'%,keyboard
            N(end+1:fianti(end))=0;
        end
        %      N=0;
        if N(1)<0
            N=-N;
        end
    else
        roadd=dox+1;
        N=exp(-(ro_in/roadd).^4);
%         g0=0;
        g0=abs(min(matgain));
    end
%    if(max(matgain)<0)
%%     g0=0;
%     velm.Lm=mode.Lm; % prendo solo il primo modo (per ora)
%     velm.vlambda=mode.vlambda;
%     velm.NQW=mode.NQW; % prendo solo il primo modo (per ora)
%     velm.Gamma_z=mode.Gamma_z;
%     velm.nindexQW=mode.nindexQW;
%     velm.fPES=mode.fPES;
%     velm.fPdif=mode.fPdif;
%     velm.E2=mode.E2; % original    
%    end
end

%N=0;
%g0=0;
pimra=-g0/cconv_velm;
calop.pimra=pimra;
%'g0', keyboard
calop.zeta=mesh.ygrid*1e4;


 %dndT0=2.1e-4* mesh.fatt_dndT;
 dndT0=mesh.dndT;

 %dndT0=2.6e-4*mesh.fCondTerTOT;
 %dndT0=3.5e-4*mesh.fCondTerTOT;
 nmed=3.248; %media di n con |E|^2
 %Ps.dndT0=dndT0*(mode.T0-300);   % valore 1D

 T_por=flipud(Deltan-1i.*ImagDeltan);%*0;
 T_ter=dndT0*Tdd;
 Tg=2*nmed*(T_ter+T_por);
 %Tg=2*nmed*(T_ter);
% 'ver T por', keyboard
% 'ver T por', keyboard

 zTf=calop.zeta; 

 dz=mesh.dz_Term/1000; 
 Ps.itutmir=itutmir;
 
 Ps.zTemp_start=dz;
 
 zT=(max(zTf)-fliplr(zTf))+dz;
 fiT=find(diff(zT)<.9*max(diff(zT)));
 Ps.dndT0=dndT0*(mean(Tdd(fiT,end)));   % valore 1D
% 'vedi T0', keyboard
 %Ps.dndT0=dndT0*(mean(Tdd(fiT,end)));   % valore 1D
 zAbs=(max(zTf)-fliplr(zTf));
 calop.zeta=zT;
global EO ABS
 EO.Es=flipud(Es);   % cambio direzione riferimento.
 EO.zT=zT;
 %'ver EO', keyboard

% parte per assorbimento free-carrier
%' abs', keyboard
ABS.zAbs=zAbs;

rFCA=dox*2/3;
[~,indox]=min(abs(NPx-rFCA));
indox=NPxQW;

elec=reshape(mode.elecABSvelm,mesh.nny,mesh.nnx);
eleccentro=elec(:,indox);
hole=reshape(mode.holeABSvelm,mesh.nny,mesh.nnx);
holecentro=hole(:,indox);

xmo=reshape(mesh.xmol,mesh.nny,mesh.nnx);
xmol=xmo(:,1);

ABS.xmol=xmol;
ABS.eleccentro=eleccentro;
ABS.holecentro=holecentro;
% ABS.zT=zT;


 % effetti termici
 
Al=mode.ABS_Apor;

'Check proper Tvelm',% keyboard
% Tvelm=Tdd(:,end);
Tvelm=Tdd(:,indox);
ABS_Texp=mode.ABS_Texp+mode.PerCoefExT*Tvelm;

if Al>0
    N0=1; % per N=1 tutte coincidono
    InDe=1/log(Al*N0+1);
    if isfield(mode,'IOLDsw')
        if mode.IOLDsw==1
            N0=mode.ABS_Ader;
        end
        InDe=(Al*N0+1)/Al;
    end
    
    ABS.eleccentro=InDe*log(Al*ABS.eleccentro+1).*(1+Tvelm/T300).^ABS_Texp;
    ABS.holecentro=InDe*log(Al*ABS.holecentro+1).*(1+Tvelm/T300).^ABS_Texp;
else
    ABS.eleccentro=ABS.eleccentro.*(1+Tvelm/T300).^ABS_Texp;
    ABS.holecentro=ABS.holecentro.*(1+Tvelm/T300).^ABS_Texp;
end
%  'qui ABS f_callVelm', keyboard

if isfield(mode,'Fat_PerCoefTemp')
    Fat_Perd=mode.Fat_Perd0+Tvelm*mode.Fat_PerCoefTemp;
    ABSh=mode.ABSh0.*Fat_Perd;
    ABSe=mode.ABSe0.*Fat_Perd;
else
    ABSh=mode.ABSh;
    ABSe=mode.ABSe;
end
ABS.e=ABSe;
ABS.h=ABSh;
 
 %alN3=(mode.holeABS.*ABSh+mode.elecABS.*ABSe)*100;
  
  
if ifp_Opt==10
  'qui ABS f_callVelm', keyboard
end  
% 'controllo T', keyboard
%  zitr=zT;
%  fitr=1;
 %Tz=dndT.*ones(size(zT));
 

 xT=NPx(fianti);
 
%  'Tg', keyboard
 ax=dox*3;
 Ps.ax=ax;
 %ax=dox*1.3;
 


 
 if T0>0
  Ps.T=Tg;
  Ps.Tx=xT;
 % Ps.Tz=zT+dz;
%   Ps.Tz=Tz;
  Ps.zT=calop.zeta;
%   ndz=3;
%   Ps.ndz=ndz;
  Ps.zitr=calop.zeta;
%   Ps.fitr=fitr;

% figure, 
% subplot(121), plot(xT,Tdd)
% xlabel('ro')
% ylabel('Delta T')
%  %'qui Temp', keyboard
%  subplot(122), plot(zT,Tdd)
%   xlabel('z')
%   pausak

if(ifp_Opt==-10)

 figure, 
 subplot(121), plot(xT,Tdd)
 xlabel('\rho')
 ylabel('Incremento Temperatura')
  %'qui Temp', keyboard
  subplot(122), plot(zT,Tdd)
   xlabel('z')
%   title('real part')
   pausak
 figure, 
 subplot(121), plot(xT,imag(Tg)*(2*k0/100))
 xlabel('\rho')
  ylabel('Delta Eps (Imag)')
  %'qui Temp', keyboard
  subplot(122), plot(zT,imag(Tg)*(2*k0/100))
     title(' LOSSES')
   xlabel('z')

  % 'qui Temp', keyboard
end

else
   Tg=0;
%    'qui Temp', keyboard

 
 end
 
%Tg=0

calop.epsT=Tg.';

calop.ro_in=ro_in';
%'ro_in', keyboard
calop.N=N;
calop.Nref=Nref;

%calop.Nref=mode.matDn;
calop.ndz=25;


  
calop.del_n_ag=del_n_ag;
calop.del_n_agInf=del_n_agInf;

if del_n_ag==0
 ianti_gui=0;
else 
 Dlai=Dlam;
 %Dlai(1:2)=Dlai(1:2)-1.3*calop.del_n_ag;
 calop.Dlam     =Dlai;
end
calop.ianti_gui=ianti_gui;

if isfield(VelmOptions,'isoga')
 Ps.isoga=VelmOptions.isoga;
else
 Ps.isoga=1;
end

%'cont matDn', keyboard
 % da riempire per SPP
 Par=0;
 Ps.Par=Par;
 calop.Ps=Ps;


%calop.N=0;
%calop.Nref=0;
%calop.epsT=0;
%calop.ianti_gui=0;
%calop.pimra=0;
%calop.del_n_ag=0;
%' ho messo tutto a ZERO in f_CallVELM',
% calop.epsT=calop.epsT*10;

if ifp_Opt==-10
  ' prima di calop in set_opt ', keyboard
end  

[St_wave,camzT,cam2_0,gain_0,lamod_0,modc,xro,fian,fPdif_0,fPES_0,fPGA,...
tyPmod_0,omP0_0,eta_eff_0,Tper,Pa,Ppol,Plot,Eo_x,Eo_y]=...
 caloptdd(MODC0,ro_campo,ro_in,fil_str,iLP,ifp_Opt,calop);

 %figure,plot(xro,cam2_0),ylabel('campo al quadrato'), keyboard
%  ' DOPO di calop in set_opt ', keyboard
 
if fianti(end)>NPxQW || ired_int==1
    'ioldE2fit'%,keyboard
    ioldE2fit=0 % 1: OLD E2 fit
else
    ioldE2fit=1 % 1: OLD E2 fit
end

% ioldE2fit=0

if ioldE2fit==1
    xro=ro_Fie/1e4;
    
    Fdox=1.3;
    % comment for 1D simulation
    fifit1=find(ro_campo>dox*Fdox);
    fifi(1)=fifit1(1);
    fifit2=find(ro_campo>dox*Fdox*1.1);
    fifi(2)=fifit2(1);
    %
    fi1=find(ro_Fie<=dox*Fdox); % spline in a region slightly larger than OX
    fi2=find(ro_Fie>dox*Fdox);
    
    % tolgo fit coda
    % fi1=1:length(ro_Fie);
    % fi2=[];
    fifit1=1:length(ro_campoF);
    
    for indMode=1:length(gain_0)
        E2tmp=cam2_0(indMode,:);    % Extracted by VELM on its mesh
        colog=polyfit(ro_campo(fifi),log10(E2tmp(fifi)),1);
        %         Etmp=E2tmp(fifi(1):end);
        %     xtmp=ro_campo(fifi(1):end);
        %     dE=diff(Etmp);
        %     dE0=dE(1:end-1).*dE(2:end);
        %     fin=find(dE0<0);
        %     Ema=Etmp(fin);
        %     xma=xtmp(fin);
        %     if Ema(1)>Ema(2)
        %         Emafin=Ema(1:2:end);
        %         xmafin=xma(1:2:end);
        %     else
        %         Emafin=Ema(2:2:end);
        %         xmafin=xma(2:2:end);
        %     end
        %     colog=polyfit(xmafin,log10(Emafin),1); % comment for 1D simulation
        
        % spline in a region slightly larger than OX
        Efit1=spline(ro_campo(fifit1),E2tmp(fifit1),ro_Fie(fi1));
        
        % Exponential fit on a region just after the ox, extended to whole
        % radius
        if length(fi2)>0
            Efit2=10.^(polyval(colog,ro_Fie(fi2)));
        else
            Efit2=[];
        end
        E2t=[Efit1 Efit2];
        Norm2=trapz(xro,2*pi*xro.*E2t);
        mode.E2(indMode,:)=E2t/Norm2;
        ca2(indMode,:)=E2tmp/Norm2;
    end
    
    
elseif ired_int==1
    %% Brings the fields to zero after a certain radius beyond
    % the oxide
    xro=ro_Fie/1e4;
%     Fdox=2;
%     Fat_ox=2;    % Same defined in ired_int
    fi1=find(ro_Fie<=dox*Fat_ox);
    fi2=find(ro_Fie>dox*Fat_ox);

    for indMode=1:length(gain_0)
        E2tmp=cam2_0(indMode,:);
        Efit1=spline(ro_campo,E2tmp,ro_Fie(fi1));
        Efit2=zeros(size(ro_Fie(fi2)));
        
        E2t=[Efit1 Efit2];
        Norm2=trapz(xro,2*pi*xro.*E2t);
        mode.E2(indMode,:)=E2t/Norm2;
        ca2(indMode,:)=E2tmp/Norm2;
    end
else
    %% Used when fianti=indexes of DD, including the PASSIVATION!
    % Be careful that N and Nref that come from matgain has set to zero
    % from the mesa to the end of passivation, otherwise dimension mismatch
    % arise in vel_DD!
    xro=ro_Fie/1e4;
    for indMode=1:length(gain_0)
        E2tmp=cam2_0(indMode,:);
        
        E2t=spline(ro_campo,E2tmp,ro_Fie);
        Norm2=trapz(xro,2*pi*xro.*E2t);
        mode.E2(indMode,:)=E2t/Norm2;
        ca2(indMode,:)=E2tmp/Norm2;
    end
    
end
 
iverfitCa=0;
if iverfitCa==1
 figure, semilogy(ro_Fie,mode.E2)
 hold on,
 semilogy(ro_campo,ca2,'.')
 semilogy(ro_campo(fifi),ca2(:,fifi),'o')
 
 pausak
 figure, plot(ro_Fie,mode.E2)
 hold on,
 plot(ro_campo,ca2,'.')
 pausak 
'campi', keyboard
end 


%for indMode=1:length(gain_0)
%    E2tmp=cam2_0(indMode,:);
%    E2tmp=E2tmp/max(E2tmp);
%    fi=find(E2tmp(fiox)<.02);
%    E2tmp(fiox(fi(1)):end)=1e-5;
%    Norm2=trapz(xro,2*pi*xro.*E2tmp);
%    mode.E2(indMode,:)=E2tmp/Norm2;
%end
%'qui dopo VELM',keyboard
cazfit=spline(camzT(:,1),(camzT(:,2)),zT);
cazfit=cazfit/cazfit(1);
finv=find(zT<min(camzT(:,1)) | zT>max(camzT(:,1)));
if length(finv)>0
 cazfit(finv)=0;
end
Inviluppo_SW=fliplr(cazfit);
velm.Inviluppo_SW=Inviluppo_SW;
St_wave(:,2)=St_wave(:,2)/max(abs(St_wave(:,2)));
velm.SW=St_wave;
%figure, semilogy(zT,cazfit,camzT(:,1),camzT(:,2))
% 'ver  cam z', keyboard, keyboard
%mode.camz=camz;


[xro2,ia,ic]=unique(NPx);

dbclear if error

%-- VELM returns mode.vlambda, mode.vGth, ...
%velm.Lm=gain_0*Pa.confztot(1);
velm.Lm=Pa.losm;
velm.dlam=Pa.delf*lambda*1000;
velm.vlambda=lamod_0*1000;
velm.Gamma_zTOT=Pa.confztot(:,2);
velm.Gamma_z=Pa.confztot(:,3:end);

velm.NQW=Pa.confztot(:,2)/Pa.confztot(:,1);
velm.nindexQW=Pa.rr;
velm.nindexRef=Pa.rr;
velm.fPES=[fPES_0];
velm.fPdif=[fPdif_0];
%' fine VELM', keyboard, keyboard

velm.E2=mode.E2(:,ic);
velm.xro=xro;

if abs(size(velm.Gamma_z,1)-mode.nmodes)>0
    fprintf('\n velm.Gamma_z is shorter than mode.nmodes!!!!\n\n')
    velm.Gamma_z(length(Pa.losm)+1:mode.nmodes,:)=velm.Gamma_z(end,:);
    velm.Lm(length(Pa.losm)+1:mode.nmodes)=velm.Lm(end);
    velm.dlam(length(Pa.losm)+1:mode.nmodes)=velm.dlam(end);
    velm.vlambda(length(Pa.losm)+1:mode.nmodes)=velm.vlambda(end);
    velm.fPdif(length(Pa.losm)+1:mode.nmodes)=velm.fPdif(end);
    velm.fPES(length(Pa.losm)+1:mode.nmodes)=velm.fPES(end);
    velm.E2(length(Pa.losm)+1:mode.nmodes,:)=velm.E2(end,:);
%     keyboard
end

Ef=Plot.Ef{1};
X=Plot.X{1};
xxf=X(:,1);

If=(Ef/max(max(Ef))).^2;
e2=exp(-2);
[du,fiL]=min(abs(mean(If,2)-e2));
ffL=xxf(fiL);
velm.If=If;
velm.theta=xxf;
velm.ffL=ffL;
velm.x=Plot.XP{1};
velm.Eoutx=Plot.E2xo{1};
if VelmOptions.ivett==1
velm.Eouty=Plot.E2yo{1};
end


%' fine VELM', keyboard, keyboard
if ifp_Opt==-10
  ' Dopo calop in f_callVelm ', keyboard
end  
warning on
return