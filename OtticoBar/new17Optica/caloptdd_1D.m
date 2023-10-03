
%
%
% function [cam2,gain,lamod,modcu,xro,fian,fPdif,fPES,fPGA,tyPmod,omP0,Ppol,Plot,Pa]=...
%           calopt(modc0,idyn,i2Ddyn,NPx,NPfi,fil_str,iLP,ifp,calopt_in)
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Input parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%     OBBLIGATORI
%  idyn   :  1, par. per dinamica, 0, solo modi
%  i2Dyn  :  1, per dinamica 2D; 0, din. 1D
%  NPx    :  # punti radiali
%  fil_str:  file struttura
%  iLP    :  =0, sol. vettoriale, =1, sol. scalare
%  ifp:      flag di controllo e visualizzazione.
%           -4(no stop), -10(1 stop), 0..3(stop crescenti)
%
%  calopt_in: struttura di 11 parametri facoltativi
%
%  1 Dlam_mod:    intervallo ricerca (default 1 nm)
%  2 mod_max:     n. max modi nell'int. di ricerca (default=0: tutti)
%  3 ro_in  :     coordinate di N (portatori QW) e epsT(lente termica)
%  4 zeta:        coordinate di N (portatori QW) e epsT(lente termica)
%  5 epsT :       lente termica e portatori (default: assenti)
%  6 N :          lente termica e portatori (default: assenti)
%  7 numodiacc:   numero modi azimutali accoppiati (default a seconda della struttura)
%  8 alim:        intervallo k trasversale modi (default a seconda della struttura)
%  9 dk1:         spaziatura modi trasversali (default a seconda della struttura)
% 10 ch_pset:     parametri da cambiare
% 11 parset :     valori parametri
% 12 Ev_or_Od :   sol. Even oppure Odd
% 13 ipolar:      soluzione cos (1) o sin(-1) o entrambe (2)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Output parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  cam2     : modulo quadro dei campi
%  xro      : coordinata radiale
%  fPdif    : costante per il buco spaziale      (eq. N)       Bs=((fPdif.*P)'*Fintv)';
%  fPES     : costante per em. spontanea nel modo (eq. P)
%  fPGA     : fattore di confinamento longitudinale
%  omP0     : perdite
%  Plot     : stuttura con varie uscite
%  Ppol     : stuttura con varie uscite

function [St_wave,camzT,cam2,gain,lamod,modcu,xro,fian,fPdif,fPES,fPGA,tyPmod,omP0,eta_eff,Tper,...
    Pa,Ppol,Plot,Eo_x,Eo_y,ord,velm1D]=caloptdd(modc0,ro_campo,NPx,fil_str,iLP,ifp,calopt_ins,iraf_int)


%whos St_wave camzT cam2 gain lamod
%whos modcu xro fian fPdif
%whos fPES fPGA tyPmod
%whos omP0 eta_eff Tper
%whos    Pa Ppol Plot
%whos Eo_x Eo_y ord velm1D

global ianti_gui iba



velm1D=0;


idyn=1;
i2Ddyn=0;
NPfi=15;

MSG = nargchk(5,16,nargin);
if length(MSG)>0
    disp(' calopt: too few parameter ')
    keyboard
end
%'calopt', keyboard
global del_n_ag pimra NPlam dlam
calopt_in=struct2cell(calopt_ins);
[Dlam,mod_max,ro_in,zeta,epsT,N,Nref,numodiacc,alim,dk1,par_in,Ev_or_Od,ipolar,...
    pasnui,iridi,ispeedi,iPERD_vol,iplan,del_n_ag,pimra,ianti_gui,box,NPlam,dlam,Ps,ndz]=deal(calopt_in{:});

%'calopt', keyboard
global inuG inuT
if length(box)>0
    rh=box;
    global PGain PcGain PRef PcRef dn_ag gprof
    global Pterm PcTerm
    PcGain=rh.PcGain;
    if inuG==0
        PGain=rh.PGain;
        PRef=rh.PRef;
        PcRef=rh.PcRef;
    end
    PcTerm=rh.PcTerm;
    if inuT==0
        Pterm=rh.Pterm;
    end
    dn_ag=rh.dn;
    gprof=rh.Gain;
end
% ' Assegno Pterm in calopt ', keyboard
%' assegno dopo', keyboard
%
% set defaults
%

if length(Dlam)==0
    iDla=1;
    Dlam_mod(1)=0;
    Dlam_mod(2)=1;
    Dlam_mod(3)=3;
elseif length(Dlam)==1
    iDla=1;
    Dlam_mod(1)=0;
    Dlam_mod(2)=Dlam;
    npuf=floor(Dlam/3);
    if npuf<3
        npuf=3;
    end
    if npuf>10
        npuf=10;
    end
    Dlam_mod(3)=npuf;
elseif length(Dlam)==2
    Dlam_mod(1)=0;
    Dlam_mod(2)=Dlam(1);
    Dlam_mod(3)=Dlam(2);
elseif length(Dlam)>=3
    Dlam_mod=Dlam;
end

%'Dla', keyboard

if length(mod_max)==0
    mod_max=100;
end

if length(epsT)~=0
    if length(ro_in)==0
        disp(' mancano i valori di ro per epsT ')
        keyboard
    end
    if length(zeta)==0
        disp(' mancano i valori di zeta per epsT ')
        keyboard
    end
end

if length(N)~=0
    if length(ro_in)==0
        disp(' mancano i valori di ro per N ')
        keyboard
    end
end

if length(ro_in)==0
    ro_in=0;
end
if length(zeta)==0
    zeta=0;
end
if length(epsT)==0
    epsT=0;
end
if length(N)==0
    N=0;
end
T=epsT;
if length(find(T~=0))==0
    %' calopt T', keyboard
    T=0;
end



inum=0;
if length(numodiacc)==0
    numodiacc=0;
    inum=1;
end

iali=0;
if length(alim)==0
    alim=0.2;
    iali=1;
end

idk1=0;
if length(dk1)==0
    dk1=0.005;
    idk1=1;
end



% verifica defaults

[z,ro,th,Sig,NQW,lambdaNot,d,radii,shav,flsv,iauto]=Lay_ct(fil_str);

% set modi accoppiati a seconda della struttura

if inum==1
    nuvet(2)=3;
    nuvet(3)=1;
    nuvet(4)=1;
    nuvet(5)=7;
    nuvet(6)=7;
    for k=[3 4 2 5 6]
        istru=find(shav==k);
        if length(istru)>0
            numodiacc=nuvet(k);
        end
    end
    istru=find(shav==3 | shav==4);
    if length(istru)>0
        ax=radii.a;
        ay=radii.b;
        ix=find(ax>0);
        iy=find(ay>0);
        lix=length(ix);
        liy=length(iy);
        if lix>liy
            rax=ax(iy)./ay(iy);
        else
            rax=ax(ix)./ay(ix);
        end
        fax=find(rax<.5 | rax>2);
        if length(fax)>0
            numodiacc=3;
        end
    end
end


% set discretizzatione modi a seconda della struttura
if iali==1
    ax=radii.a;
    ay=radii.b;
    ix=find(ax>0);
    iy=find(ay>0);
    miax=min(ax(ix));
    miay=min(ay(iy));
    mis=min([miax miay]);
    alv=[.05,.1,.15,.2,.28];
    rav=[10,5,3,1.5,.5];
    duf=find((mis-rav)>0);
    fia=duf(1);
    alim=alv(fia);
end

if length(alim)==3
    alim_in(1)=alim(1);
    alim_in(2)=alim(2);
    alim_in(3)=alim(3);
elseif length(alim)==1
    alim_in(1)=alim;
end



% set discretizzatione modi a seconda della struttura
if idk1==1
    ax=radii.a;
    ay=radii.b;
    miax=max(ax);
    miay=max(ay);
    mas=max([miax miay]);
    dkv=[.001,.003,.005,.007,.01];
    rav=[10,5,3,1.5,.5];
    duf=find((mas-rav)>0);
    fia=duf(1);
    dk1=dkv(fia);
end


%%%%%%%%%%
%
%disp(' debug provvisorio in calolpt')
%alim=.25;
%dk1=.015;




i2D=3;
iffcut=0;


iany=1;  % 1: anisotropia planare, 2: confinata (effetto elettro-ottico )
if isfield(Ps,'iany')
    iany=Ps.iany;
end
%'iany', keyboard
%'cambiato iany !!!!!'
%iany=0  % 1: anisotropia planare, 2: confinata (effetto elettro-ottico )
ianys=0;   % strain


if iLP==1
    ianys=0;
    iany=0;
end

global fattany
fattany=-1;
fattany=1;

%' in calopt ', keyboard
iraff=0;


nmasce=mod_max;


if length(Ev_or_Od)==0
    if mod_max==1
        Ev_or_Od='Even';
    elseif mod_max>1
        Ev_or_Od='Both';
    end
end


if length(ipolar)==0
    ipolar=2;  % -1, 1, 2 (entrambe), 0 (entrambe accoppiate)
    if iLP==1
        istru=find(shav>0);
        inocir=1;
        if length(istru)>0
            inocir=length(find(shav(istru)~=1));
        end
        if inocir==0
            ipolar=1;
        end
    end
end

%keyboard
%if exist('cce')==0
% cce=0;
% save cce cce
%end
%'cal', keyboard


global pasnu

if length(pasnui)==0
    if iLP==1
        if imag(cce)==0
            ipolar=1;
            pasnu=1;
        else
            ipolar=0;
            pasnu=1;
        end
        
        if abs(cce)==0
            ipolar=1;
            pasnu=2;
        end
    end
else
    pasnu=pasnui;
end


global ispeed irid_bas

if length(iridi)==0
    irid_bas=1;
else
    irid_bas=iridi;
end


if length(ispeedi)==0
    if ipolar~=0
        ispeed=1;
    else
        ispeed=0;
    end
else
    ispeed=ispeedi;
end

%'ca ri',keyboard

fileeps='Gres_jo2';

xroI=ro_campo;
%xroI=NPx;
r_pil=[];



global igamveb igamveu GGbext GGuext
global ilo

igamveb=1;   %=1, gamma vero,  =0, riflessione GGext
igamveu=1;
GGbext=1;
GGuext=1;


npfi0=NPfi;

if npfi0==25
    
    lfi_inp(1)=4*npfi0;
    lfi_inp(2)=200;
    
else
    
    lfi_inp(1)=4*npfi0;
    lfi_inp(2)=2*lfi_inp(1);
    %lfi_inp(2)=400;
    
end
%lfi_inp(1)=4*npfi0;
%lfi_inp(2)=200;


fimaxi=pi;

fimaxi=2*pi;
%' ianti calop', keyboard
%iplan=1;
%if ianti_gui==1
% iplan=0;
%end


if length(alim)==1
    nK_dis(1)=ceil(alim/dk1(1));
elseif length(alim)==3
    nK_dis(1)=ceil(alim(2)/dk1(1));
    nK_dis(2)=ceil((alim(3)-alim(2))/dk1(2));
end

if ifp==10
    par_in
    pausak
    
else
    
    
    %   disp(' call calopt dentro '),
    %   keyboard, keyboard
    
    
    t0=clock;
    
    if exist('iraf_int')==0
        iraf_int=0;
    else
        if length(iraf_int)==0
            iraf_int=0;
        end
    end
    global NMO_in
    NMO=NMO_in;
    
    
    
    Dla=-.005;
    nomes='dp';
    ivfre0=0;
    veac=1;
    r_pil=par_in{1};
    global Ps
    global Lmax irec_fie iord_long num_long
    Lmax=100000;
    iord_long=Ps.iord_long;
    
    if length(T)>1
        zOpt=sum(flsv([1:end-1],2).*d(1:end-1));
        fi=find(zeta<zOpt);
        Topt=T(:,fi);
        zopt=zeta(fi);
    else
        Topt=T;
        zopt=zeta;
    end
    %Topt=T;
    %zopt=zeta;
    % ' VEL DD', keyboard
    [St_wave,camzT,Eqw,Eout,xro,fian,lambda,delf,gain,ord,nrAzim,Cu,Pa,ADom,...
        ADcm,Plot,Ppol,Dla_out,imod_err,velm1D]=...
        vel_DD_1D(Nref,N,Topt,ro_in,zopt,xroI,fimaxi,lfi_inp,par_in,r_pil,fil_str,...
        ipolar,Ev_or_Od,nmasce,numodiacc,Dlam_mod,nK_dis,alim,...
        iplan,iraff,idyn,i2D,iLP,iany,ianys,ifp,Dla,fileeps,nomes,ivfre0,1,ndz,veac);
    %   ' dopo VEL', keyboard
    
    if isfield(Ps,'i1D')
        i1D=Ps.i1D;
    else
        i1D=0;
    end
    
    if i1D==1
        lamod=0;
        gain=0;
        modcu=0;
        cam2=0;
        % camzT=0;
        fPdif=0;
        fPES=0;
        fPGA=0;
        tyPmod=0;
        omP0=0;
        eta_eff=0;
        Tper=0;
        Eo_x=0;
        Eo_y=0;
        %'variabilie da vedere', keyboard
        return
    end
    
    
    Pa.delf=delf;
    
    tim=etime(clock,t0);
    ore=fix(tim/3600);
    minuti=fix((tim-ore*3600)/60);
    secondi=(tim-ore*3600-minuti*60);
    format short
    if ifp==-10
        disp(' Partial elapsed time (hours, minutes, seconds) = ');
        disp([ ore minuti secondi]);
    end
end


if idyn==1
    sE=size(Eqw);
    if i2D==0
        modcv=[1:sE(2)];
    else
        if length(sE)==3
            modcv=[1:sE(3)];
        else
            modcv=1;
        end
    end
    %   modc0=[1:sE(end)];
    modc=modcv;
    if modc0~=0
        if length(modcv)>length(modc0)
            modc=modc0;
        end
    end
    %    modcu=modcv;
    modcu=modc;
    
    %' modc ', keyboard
    
    
    
    %   Gam_u =Pa.Gam_u;
    %   Wvmo =Pa.W;
    %   fa4mo =Pa.famod;
    %   czmo =Pa.cz;
    %   fqwmo=Pa.fqw;
    %   Lfmo =Pa.Lf;
    %   Tfmo =Pa.Tf;
    %   Tfimo=Pa.Tfinf;
    %   gtmo =Pa.gt;
    %   fatqw=Pa.fatqw;
    %   d_at =Pa.dat;
    %   NQW  =Pa.NQW;
    %   rr   =Pa.rr;
    %   avero=Pa.avero;
    
    Tefv=Pa.trasm;
    
    
    rr   =Pa.rr;
    rvg   =Pa.rmed;
    %   fqwmo=Pa.fqw;
    d_at =Pa.dat*1e-6;    % in metri
    NQW=Pa.NQW;
    dat_tot=d_at*NQW;
    tyPmodu=Pa.type;
    %   'calopt rr', keyboard
    
    tyPmod=tyPmodu(modc);
    
    taut=Pa.taut;
    tauu=Pa.tauu;
    taub=Pa.taub;
    %   alca=Pa.alca;
    pvol=Pa.pvol;
    
    confztot=Pa.confztot(2);
    confz=Pa.confztot(1);
    
    %   losm=Pa.losm/confz*confztot/NQW;
    losm=Pa.losm;
    
    
    %   fatqw=mean(fqwmo);
    %   losm=;
    %   Lefv=Lfmo;
    %   Tefv=gtmo;
    %   Tefiv=Tfimo;
    
    %   confz=czmo;
    %   confztot=czmo.*NQW.*fqwmo;
    frism=delf(modc);
    %global ijos tauN0 hbar c htom
    global ijos  hbar c_light htom tauN0
    
    
    om0=2*pi*1e6*c_light/lambda;
    htom=hbar*om0;
    lamod=lambda./(1+frism);
    
    %' calopt dope ', keyboard
    
    % set discretizzazione spaziale
    
    differen
    
    % end set discretizzazione spaziale
    %' % end set discretizzazione spaziale '
    
    %'cam2 ', keyboard
    
    
    
    
    if i2Ddyn==1
        
        for k=modc
            F=abs(reshape(Eqw(ifidif,ifisho,k),length(ifidif),lfs)).^2;
            F1=abs(reshape(Eqw(:,ifisho,k),length(xro),lfs)).^2;
            Sint=spacintv(F1,xdx',lfp,pesf,Nsett);
            cam2(k,:)=reshMv(F)'/Sint;
            Fx=abs(reshape(Eout.x(ifidif,ifisho,k),length(ifidif),lfs)).^2;
            Fy=abs(reshape(Eout.y(ifidif,ifisho,k),length(ifidif),lfs)).^2;
            Sint=spacintv(Fx+Fy,xdx1',lfp,pesf,Nsett);
            Fx=reshape(Eout.x(ifidif,ifisho,k),length(ifidif),lfs);
            Fy=reshape(Eout.y(ifidif,ifisho,k),length(ifidif),lfs);
            Eo_x(:,:,k)=Fx/sqrt(Sint);
            Eo_y(:,:,k)=Fy/sqrt(Sint);
        end
        
    else
        Eo_x=0;
        Eo_y=0;
        if i2D==3
            
            for k=modc
                F=abs(reshape(Eqw(ifidif,ifisho,k),length(ifidif),lfs)).^2;
                F1=abs(reshape(Eqw(:,ifisho,k),length(xro),lfs)).^2;
                Sint1=spacintv(F1,xdx',lfp,pesf,Nsett);
                Sint=spacintv(F,xdx1',lfp,pesf,Nsett);
                cdp=pesf*F';
                cdp1=pesf*F1';
                cam2(k,:)=cdp/Sint1*Nsett/(2*pi);
                cam2v(k,:)=cdp1/Sint*Nsett/(2*pi);
            end
            
        else
            
            F=abs(Eqw(ifidif,:)).^2;
            F1=abs(Eqw).^2;
            Sint=diag(1./(2*pi*xdx*F));
            Sint1=diag(1./(2*pi*xdx1*F1));
            cam2=Sint1*F';
            cam2v=Sint*F1';
            
        end
    end
    
    ifi=0;
    omP0=losm(modc)';
    fPGA=confztot;
    
    %      'qui fPGA', keyboard
    %      Ttot=Tefv(modc)';  % per uscita solo dall'alto
    %      Nu_tot=fPGA.*tauu(modc)'/dat_tot;  % per uscita solo dall'alto
    %      fPdif_old=tauN0*4*1e9./(c/rr*hbar*om0.*(1+frism).*Ttot)*1e-24;
    %      fPdif=tauN0*1e9./(hbar*om0.*(1+frism))*1e-24.*Nu_tot;
    Lef=Pa.Lf;
    CP=2*Lef/(c_light/rvg);
    Ttc=CP./tauu;
    %      Pp1=(tauu./Lef)';
    %      fPdif=tauN0*2*Pp1(modc)*1e9./(hbar*om0.*(1+frism))*1e-24;
    %      fPdif=tauN0*Pp1(modc)*1e9./(hbar*om0.*(1+frism))*1e-24;
    
    fatqw=Pa.fatqw;
    vr=c_light/rvg*1e2;
    
    Teff=Pa.Teff;
    Ppn=1./(Teff'*vr)*fatqw;
    
    vr1=c_light/rr*1e2;
    %      vr1=c_light*1e2;
    %      'modifica fPdif !!! '
    
    %      Ppn=1./(Teff'*vr1);
    
    %      Ppn=1./(Teff'*vr1)*fatqw*NQW;
    %      Ppn=rr./(Teff'*vr1)*fatqw; % RIFLETTERE ANCORA SU QUESTO!!! fatqw!!
    
    %      Ppn=rr./(Teff'*vr1); % RIFLETTERE ANCORA SU QUESTO!!! fatqw!!
    
    Ppnli=1./(Teff)';
    %       fPdif=tauN0*Ppn(modc)*1e9./(hbar*om0.*(1+frism))*1e-24/1000;
    
    fPdif=1./(Teff.*hbar*om0.*(1+frism')); %-- tibaldi
    global fPnli
    %     'fPdif', keyboard
    
    %      Teq=1./(1./Tefv+1./Tefiv);
    %      taP=2*Lefv/(c/rr.*Teq).*fatqw;
    %      taP=2*Lefv/(c/rr.*Teq).*fatqw;
    
    taP=tauu(modc)';
    fPnli=1e-3*Ppnli./(hbar*om0.*(1+frism)).*taP;
    'fPdif', keyboard
    fPES=hbar*om0*(1+frism)./taP;
    %      fPES=tauN0*hbar*om0*(1+frism)./taP*1e3;
    %      fPES=tauN0*hbar*om0*(1+frism)./taP*1e3.*fatqw;
    eta_eff=(taut./tauu)';
    Tper{1}=taut';
    Tper{2}=tauu';
    Tper{3}=taub';
    Tper{4}=pvol';
    %      Tper{5}=alca';
    %' qui pvol calotp taut ', keyboard
    %' qui pvol calotp taut ', keyboard
    %' qui pvol calotp taut ', keyboard
    
    if ifp>=0
        ' %% fine set costanti',
        keyboard
    end
    
    % function [cam2,gain,lamod,modcu,xro,fian,fPdif,fPES,fPGA,tyPmod,omP0,...
    % eta_eff,Tper,Ppol,Plot,Pa,Eo_x,Eo_y]=calopt
    
    
    
else  % idyn
    
    [cam2,xro,gain,lamod,modc,fPdif,fPES,fPGA,omP0]=deal(0);
    
end

'fine calopt dd', keyboard
