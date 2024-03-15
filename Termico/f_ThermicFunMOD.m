function [Tdd,T3Dt,PTherm,T_Contributi,fattore_correttivo]=f_ThermicFunMOD(Tprec,mesh,mode,StrTT,IPLOT)

%clear all
%close all
%Tprec=0;
if nargin==4
    IPLOT=2;
end
ifig=abs(mode.iTfig);     % =1 plot risultati intermedi
% ifig=1;     % =1 plot risultati intermedi
%load Dati_termici_esempio

% load output_HalfVCSEL_p_cy

%'cont Ter', keyboard
%   function [Tdd,T3D]=f_ThermicSimulator(Tprec,mesh,mode)

Temp0=mode.T300;
T_Peltier=mode.T0;
Exp_Temp0=mode.Exp_Temp0;

%  'Tempa', keyboard
%  Tprec=T3D;

if isfield(mode,'ioldTemp')==1
    iold=mode.ioldTemp;
else
    iold=1;
end


if iold==0
    Ti=((T_Peltier+Tprec)/Temp0).^Exp_Temp0;
elseif iold==2
    Et0=-1.3
    Et=Et0+2e-3*(Tprec+T_Peltier-Temp0);
    %Et=-1.3+4e-3*Tprec;
    Ti=((T_Peltier+Tprec)/Temp0).^Et;
elseif iold==1
    Ti=(1+Tprec/Temp0).^Exp_Temp0;
elseif iold==-1
    Ti=1;
end


%Ti=1;
Tbu_dd=StrTT.Tbuf_dd;
Z_sub=StrTT.Tbuf_dd;
% 'split substrate', keyboard
Tbuf=StrTT.Tbuf;
if abs(Tbuf-Tbu_dd)>.2
    DelZet=Tbuf-Z_sub;
    StrTT.Tdbr_inf= StrTT.Tdbr_inf-DelZet;
    Z_sub=Tbuf;
end




Tbu_dd=0;

%-- inizio Tibaldi
rho=mesh.xgrid;
z=mesh.ygrid(2:end-1);   % per eliminare i contatti
%'sim', keyboard
rho=rho*1e4;
zf=z*1e4-Tbu_dd;

Jouled=reshape(mode.HeatJoule,mesh.nny,mesh.nnx);
Joule=Jouled(2:end-1,:);
%ind=find(zf>7);
%Joule(ind,:)=0;

Joule=Joule.';



Tho=reshape(mode.HeatThomson,mesh.nny,mesh.nnx);
Tho=Tho(2:end-1,:);
Tho=Tho.';


Rec_nrd=reshape(mode.HeatRec_nr,mesh.nny,mesh.nnx);
Rec_nr=Rec_nrd(2:end-1,:);
Rec_nr=Rec_nr.';

Rec_nrd=reshape(mode.HeatRec_Cap,mesh.nny,mesh.nnx);
Rec_Cap=Rec_nrd(2:end-1,:);
Rec_Cap=Rec_Cap.';

Rec_nrd=reshape(mode.HeatRec_RAD,mesh.nny,mesh.nnx);
Rec_RAD=Rec_nrd(2:end-1,:);
Rec_RAD=Rec_RAD.';

Rec_nrd=reshape(mode.HeatRec_13,mesh.nny,mesh.nnx);
Rec_srhAu=Rec_nrd(2:end-1,:);
Rec_srhAu=Rec_srhAu.';

% 'Rec', keyboard

if(not(mode.Oflg))
    mode.HeatOptAbs=zeros(size(mode.HeatJoule));
end

OptAbsd=reshape(mode.HeatOptAbs,mesh.nny,mesh.nnx);
OptAbs= OptAbsd(2:end-1,:);
OptAbs=OptAbs.';

%'Thomson',  keyboard

%-- fine Tibaldi


T_DD=max(zf);




ro_pass=StrTT.ro_pass; % questi ro sono tutti radiali (um)
ro_met=StrTT.ro_met;
ro_mesa=StrTT.ro_mesa;
ro_max=StrTT.ro_max;
Rox=StrTT.Rox;

ipass=1;
if ro_pass==ro_mesa
    ipass=0;
end




% spessori z

%Tbuf=5;

fsplit=.8;	% punto di divisione nel buffer per ridurre funzioni di base

io=0;
if io==1
    Tbu_dd=0;
    Tbuf=350; % spessore buffer (um)
    TBuf_DD=1;
    Tdbr_inf=4.7; % spessore specchio dbr inferiore (um)
    Tdbr_sup=2.721;  % spessore specchio dbr superiore

    Tcentrale=T_DD-Tdbr_inf-Tmetallo-Tdbr_sup-TBuf_DD;
    Tcav=Tcentrale;
end

Tmetallo=.2; % spessore contatto metallico

Tbuf=StrTT.Tbuf;
Tdbr_inf=StrTT.Tdbr_inf;
Tdbr_sup=StrTT.Tdbr_sup;
Tdbr_sup=StrTT.Tdbr_sup;
Tcav=StrTT.Tcav;

dZ_ther_dd=Tbuf-Z_sub;

%'spessori Temp', keyboard

% parametri sorgente di calore qtot
Ri0=10; % semiasse radiale dell'ellisse
El0=0.5; % distorsione dell'ellisse
icir=1; % se 1 e` circolare
NPfi=20; %numero punti azimutali (caso non circolare); si usano 4*NPfi+1

%icir=0;


iads=1; % aggiunge strato relativo a isplit in z
%imedia=input(' Media cond. termica? [0/1] ' )
imedia=0;  % fa media conducibilita termica
itras=0;  %media trasversale
iaria=0;   % aggiunge strato d'aria sopra la struttura
th_aria=10;  % usato solo se iaria =1



iKcost=0; % se 1 considera costante a tratti la conducibilita`, se no mette l'autoconsistenza con T
iploCond=1;
if ifig==0
    iploCond=0;
end

tic



nMax=2; % armoniche azimutali

nu_num_max=nMax;  % caso circolare
if icir==1
    nu_num_max=1;  %1 caso circolare
end



ispliro=1;  % divide intervallo in rho per migliorare set funzioni


roh=[ro_met ro_mesa ro_pass ro_max];   % raggi interfacce
if ipass==0
    roh=[ro_met ro_mesa ro_max];   % raggi interfacce
end


if ispliro==1
    IPspli=2;   % intervallo in rho da dividere
    % roh=sort([roh roh(IPspli)*[.8 2]]);
    roh=sort([roh roh(IPspli)*[ 2]]);
end

Tvet=[Tdbr_inf Tcav Tdbr_sup];

imet=1;  % includo o meno contatto nel problema termico
if Tmetallo==0
    imet=0;
end




if imet==1
    Th=[Tbuf*fsplit Tbuf*(1-fsplit) Tvet Tmetallo];  %spessori in z
else
    Th=[Tbuf*fsplit Tbuf*(1-fsplit) Tvet];  %spessori in z
    Tmetallo=0;
end
iPa=3;

Lzd=sum(Th(1:end-1));

zza=[0 cumsum(Th)];
xxa=[0 roh];

z_act=sum(Th(1:iPa))+Th(iPa+1)/2;



Nr=50;
xx1=linspace(0,Rox*3/2,Nr);
xx2=linspace(Rox*3/2,roh(1),Nr);
xx3=linspace(roh(1),roh(2),Nr);
xx4=linspace(roh(2),roh(end),Nr);
xx=[xx1 xx2(2:end-1) xx3(2:end-1) xx4];

Nr=50;
xx1=linspace(0,Rox*3/2,Nr);
xx2=linspace(Rox*3/2,roh(1),Nr);
xx3=linspace(roh(1),roh(2)*1.2,Nr);
xx4=linspace(1.2*roh(2),roh(end),Nr);

xx=[xx1 xx2(2:end-1) xx3(2:end-1) xx4];



%L=60;  %numero funzioni Bessel espansione Temperature
%J=30;  % numero funzioni trigonometriche temperatura


nomet='provaAT';



% load ParT  % da settare e sistemare; per ora vecchia struttura


fian=linspace(0,2*pi,4*(NPfi)+1);

ang=fian(1:NPfi+1); % [1 * 26] ,il 26^ è pi/2


% e' settato in calopt !!!!!!
%if (i2Dyn==1 & versione_top ~=0) % ------ > 2D misto



if nu_num_max==1
    ang=0;
    NPfi=0;     %circolare
end


% e' settato in calopt !!!!!!
%if (i2Dyn==1 & versione_top ~=0) % ------ > 2D misto

nu_num=[1:nu_num_max];
nu_vett=[0 nu_num(1:end-1).*2];



if length(ang)>1
    difan=diff(ang);
else
    difan=pi/2;
end


zzd=zza(2:end-1);

if iaria==1
    zza=[zza zza(end)+th_aria];
    zzd=zza(2:end-1);
    'aria',
    %  keyboard
end


fCond=mesh.fCondTer;
fCondZ=mesh.fCondTerZ;

Bcond0=StrTT.Bcond;
Bcond.CondZc=Bcond0.CondZc*fCondZ; % cladding, from Markus, for xmol almost 0.2 (cladding + QW quindi faccio media)
Bcond.CondTc=Bcond0.CondZc*fCond; % cladding, from Markus, for xmol almost 0.2 (cladding + QW quindi faccio media)

Bcond.CondZm=Bcond0.CondZm*fCondZ; % transverse mirror, theoretical calculation Main_PLOT_ThermalCond_Vs_xmol
Bcond.CondTm=Bcond0.CondZm*fCond; % longitudinal mirrors, theoretical calculation Main_PLOT_ThermalCond_Vs_xmol

Bcond.CondZb=Bcond0.CondZb*fCond; % substrate, Markus, for xmol almost 0.0

% Bcond.Cond_air=Bcond0.Cond_air*fCondZ; % air
% Bcond.CondZmet=Bcond0.CondZmet*fCondZ; % metal
% Bcond.Condpas=StrTT.Bcond.Condpas*fCondZ; % passivation

Bcond.Cond_air=Bcond0.Cond_air; % air
Bcond.CondZmet=Bcond0.CondZmet; % metal
Bcond.Condpas=StrTT.Bcond.Condpas; % passivation

iLith=mode.flgBTJ_lithographic;

% 'qui termico', keyboard

CondZc=Bcond.CondZc;
CondTc=Bcond.CondTc;

CondZm=Bcond.CondZm;
CondTm=Bcond.CondTm;

CondZb=Bcond.CondZb;

CondZmet=Bcond.CondZmet;
Cond_air=Bcond.Cond_air;
Cond_air_side=Cond_air;
if StrTT.TrenchPass==1
    Cond_air_side=CondTm;
end

if isfield(Bcond,'Condpas')==1
    Cond_pas=Bcond.Condpas;
    %  ipass=1;
else
    Cond_pas=Cond_air;
end

if iLith==1
    Cond_air_sideZ=CondZm;
    Cond_air_sideT=CondTm;
    Cond_pasZ=CondZm;
    Cond_pasT=CondTm;
    Cond_pas=Cond_air;
else
    Cond_air_sideZ=Cond_air;
    Cond_air_sideT=Cond_air;
    Cond_pasZ=Cond_pas;
    Cond_pasT=Cond_pas;

    Cond_air_sideZ=Cond_air_side;
    Cond_air_sideT=Cond_air_side;
    %    Cond_pasZ=Cond_air_side;
    %    Cond_pasT=Cond_air_side;




end


condz1=flipud([Cond_air;    CondZm; CondZc; CondZm;  CondZb]);
condz2=flipud([CondZmet;  CondZm; CondZc; CondZm;  CondZb]);
condz3=flipud([Cond_air_sideT; Cond_air_sideZ; CondZc; CondZm; CondZb]);

condt1=flipud([Cond_air; CondTm; CondTc; CondTm;  CondZb]);
condt2=flipud([CondZmet; CondTm; CondTc; CondTm; CondZb]);
condt3=flipud([Cond_air; Cond_air_sideT; CondTc;  CondTm; CondZb]);
condt3=flipud([Cond_air_sideT; Cond_air_sideT; CondTc;  CondTm; CondZb]);


if ispliro==0
    if ipass==1
        condz23=flipud([Cond_pas; Cond_pasZ; CondZc;  CondZm; CondZb]);
        condz=[condz1 condz2 condz23 condz3];
        condt23=flipud([Cond_pas; Cond_pasT;  CondTc;  CondTm; CondZb]);
        condt=[condt1 condt2 condt23 condt3];
    else
        condz=[condz1 condz2 condz3];
        condt=[condt1 condt2 condt3];
    end
else
    condz13=flipud([CondZmet; CondZm; CondZc;  CondZm; CondZb]);
    condz23=flipud([Cond_pas; Cond_pasZ; CondZc;  CondZm; CondZb]);
    condz33=flipud([Cond_air_side; Cond_air_sideZ; CondZc;  CondZm; CondZb]);

    condt13=flipud([CondZmet; CondTm; CondTc;  CondTm; CondZb]);
    condt23=flipud([Cond_pas; Cond_pasT; CondTc;  CondTm; CondZb]);
    condt33=flipud([Cond_air_side; Cond_air_sideT; CondTc;  CondTm; CondZb]);

    if ipass==1
        condzT=[ condz23 condz33];
        condtT=[ condt23 condt33];
        condzT=[ condz23 condz33];
        condtT=[ condt23 condt33];
    else
        condtT=[condt13 condt33];
        condzT=[condz13 condz33];
        condtT=[condt33];
        condzT=[condz33];
    end
    condz=[condz1 condz2 condzT condz3];
    condt=[condt1 condt2 condtT condt3];

end

if imet==0
    condz=condz(1:end-1,:);
    condt=condt(1:end-1,:);
end

if iaria==1
    condt=[condt; ones(size(condt(end,:)))*Cond_air_side ];
    condz=[condz; ones(size(condt(end,:)))*Cond_air_side ];
    %'aria cond', keyboard
    %condz=
end

if iads==1
    condt=[ones(size(condt(end,:)))*CondZb; condt ];
    condz=[ones(size(condt(end,:)))*CondZb; condz ];
end

if itras==1
    condt=repmat(condt(:,1),size(condz(1,:)));
    condz=repmat(condz(:,1),size(condz(1,:)));
    %  condt=condt(:,1);
    %  condz=condz(:,1);
end
%  'cond', keyboard

if imedia==1
    %     condzm=mean(mean(condz));
    %     condtm=mean(mean(condt));
    %     condz(:,:)=condzm;
    %     condt(:,:)=condtm;
    if iads==1
        condzm=mean(condz(2:end,1))
        condtm=mean(condt(2:end,1))
    else
        condzm=mean(condz(:,1))
        condtm=mean(condt(:,1))
    end

    condz(:,:)=condzm;
    condt(:,:)=condtm;
    %  condz=repmat(condz(:,1),size(condz(1,:)));
    %  condt=repmat(condt(:,1),size(condz(1,:)));

end



Lvd=diff(zza);

finz=find(Lvd>0);
finza=[finz finz(end)+1];
zza=zza(finza);
Lv=Lvd(finz);
Nv=3*ones(1,length(Lv));

Nv(2:end-1)=8;
%Nv(end)=12;
Nv=5*ones(1,length(Lv));
%Nv(1)=3;

fiCav=find(Th==Tcav);
zz=[];
dz=[];
for ki=1:length(Lv)
    xini=zza(ki);
    L=Lv(ki);
    Nquad=10*Nv(ki)+1; %-- conservativo
    if ki==fiCav
        Nquad=50*Nv(ki)+1; %-- conservativo
    end
    [nodes,weights]=quadad('legen',1,Nquad);
    % mapping from [-1,+1] and [xi,xr]
    xi=L/2*nodes+L/2+xini;
    wx=L/2*weights;
    %    'ver',keyboard
    fiZ{ki}=length(zz)+[1:length(xi)];
    zz=[zz xi];
    dz=[dz wx];

end
%'cont fine', keyboard



[f,fp,fs,Nfun,MI,MIs,MIp]=SEM_Nver(Lv,Nv,zz);

%[f_int,fp_int]=SEM_Nver(Lv,Nv,zzd);


%' proma di iLO', keyboard

% S=sin(kdisL'*zz');
% Su=sin(kdisL'*(zeta'+thick_bu));
S=f;
s2=f;
s2s=fp;
L=Nfun;
I=L;

if iKcost==1
    for kj=1:L
        for ki=1:L
            SIm{ki,kj}=squeeze(MI(kj,ki,:))';
            SIm2{ki,kj}=squeeze(MIp(kj,ki,:))';
        end
    end
end

if ispliro==0
    Nvr=[10 10 5];
    Nvr=[4 6 16];
    of=[0 0 0];
else
    Nvr=[4 6 6 16 16];
    if ipass==1
        of=[0 0 0 0 0];
    else
        of=[0 0 0 0 ];
    end
    nv0=10;
    Nvr=nv0*ones(size(of));
    Nvr(end)=2;
end

% forse conviene mettere in testa il numero di funzioni che vuoi usare su ciascun patch


    


Rv=diff(xxa);

% 'Rv', keyboard

fimag=find(Rv>0);
Rv=Rv(fimag);
Nvr=Nvr(fimag);
of=of(fimag);

tic
% ioldB=input(' vecchio metodo = ');


%  Nvr=[30];
%  Rv=xxa(end);
% [fb,xb,wb,NfunB,MIB,MIpB]=SEM_BasM(Rv,Nvr);


%[f_int,fp_int]=SEM_Nver(Lv,Nv,zzd);


%' proma di iLO', keyboard

% S=sin(kdisL'*zz');
% Su=sin(kdisL'*(zeta'+thick_bu));




icnu=0;
Jsu=0;
for nua=nu_vett
    icnu=icnu+1;
    [fb,fbp,fbs,xb,wb,NfunB,MIB,MIpB,MI1B,fiX]=SEM_Bessel(Rv,Nvr,of,nua);
    fM{icnu}=fb;
    fTrr{icnu}=fbp;
    fTrf{icnu}=fb*diag(1./xb)*nua;
    J=NfunB;
    Js(icnu)=J;
    Jsu=Jsu+J;
    if iKcost==1
        for jp=1:J
            for jj=1:J
                JIm{icnu,jp,jj}=squeeze(MIB(jp,jj,:));
                JImp{icnu,jp,jj}=squeeze(MIpB(jp,jj,:));
                JIm1{icnu,jp,jj}=squeeze(MI1B(jp,jj,:));
            end
        end
    end
end

R=[];
% angi=ang+difan(1)/2;
angi=ang;
Coa=[];
for knu=1:length(nu_vett)
    Bd0=fM{knu}';
    Bd=[];
    for kfi=1:NPfi+1
        Bd=[Bd; Bd0*cos(nu_vett(knu)*angi(kfi))];
    end
    R=[R Bd];
end
Ru=R;

%  xxd=[0 xb];
%  fid=find(diff(xxd)>0);
%  xb=xb(fid);
%fid=1:length(wb);

ro_out=xb;
xx=ro_out;


xdx=wb*2*pi;
xdx0=wb;
if length(nu_vett)>1
    xdx=[wb/2 repmat(wb,1,NPfi-1) wb/2]*difan(1)*4;
end
%  'Fcotn', keyboard
% if iKcost==0
Fcondt=zeros(length(zz),length(xx));
Fcondz=Fcondt;
for kx=1:length(fiX)
    px=fiX{kx};
    %   'px',keyboard
    for kz=1:length(fiZ)
        pz=fiZ{kz};
        Fcondt(pz,px)=condt(kz,kx);
        Fcondz(pz,px)=condz(kz,kx);
    end
end

Fcondt0=Fcondt;
Fcondz0=Fcondz;

%  'Fcotn', keyboard
Fcondt=Fcondt0.*Ti';
Fcondz=Fcondz0.*Ti';

[xxc,zzc]=meshgrid(xx,zz);

if iploCond==1
    figure, surf(xxc,zzc,log10(Fcondt)), view(2), shading flat,
    title(' Transverse Thermal conducibility')
    a=axis;
    a(3)=Tbuf-2;
    a(4)=max(max(zzc))+1;
    axis(a)
    xlabel('\rho, \mum'),ylabel('z, \mum')
    colorbar
    pausak

    figure, surf(xxc,zzc,log10(Fcondz)), view(2), shading flat,
    title(' Longitudinal Thermal conducibility')
    a=axis;
    a(3)=Tbuf-2;
    a(4)=max(max(zzc))+1;
    axis(a)
    xlabel('\rho, \mum'),ylabel('z, \mum')
    colorbar
    pausak
end
% 'cond', keyboard

icnu=0;
for nua=nu_vett
    icnu=icnu+1;
    fTr=fTrr{icnu};
    fTf=fTrf{icnu};
    fT=fM{icnu};
    Mulr0=fTr*diag(xdx0);
    Mulf0=fTf*diag(xdx0);
    Mulx0=fT*diag(xdx0);
    Mulz0=f*diag(dz);
    Mulpz0=fp*diag(dz);
    %   Mulx=repmat(Mulx0,L,1);
    %   Mulz=repmat(Mulz0,Js(icnu),1);
    pJL=L*Js(icnu);
    ic=0;
    clear MLi
    %   figure
    for kx=1:Js(icnu)
        for kz=1:L
            ic=ic+1;
            MLtr=Mulz0*diag(f(kz,:))*Fcondt*diag(fTr(kx,:))*Mulr0';
            MLtf=Mulz0*diag(f(kz,:))*Fcondt*diag(fTf(kx,:))*Mulf0';
            MLz=Mulpz0*diag(fp(kz,:))*Fcondz*diag(fT(kx,:))*Mulx0';
            Msom=MLtr+MLtf+MLz;
            %     semilogy(abs(Msom)), pausak
            MLi(ic,:)=reshape(Msom,pJL,1);
        end
    end
    MLm{icnu}=MLi;

    %'quib', keyboard
end  %nua

%end %iKcost



%'fine iKost', keyboard


J=Jsu;
pJL=J*L;

icnu=0;
for nua=nu_vett
    icnu=icnu+1;
    icont=0;

    Ji=Js(icnu);

    if iKcost==1
        clear ML
        for jp=1:Ji
            %     kra=kdisJ(jp);

            for ip=1:L

                icont=icont+1;
                icont_int=0;
                clear Mdu

                for jj=1:Ji


                    JI=JIm{icnu,jp,jj};
                    JIp=JImp{icnu,jp,jj};
                    if nua~=0
                        JI1=JIm1{icnu,jp,jj}*nua^2;
                    else
                        JI1=0;
                    end
                    %         krb=kdisJ(jj);

                    for ii=1:L

                        icont_int=icont_int+1;
          	          SI=SIm{ip,ii};
        	          SI2=SIm2{ip,ii};
         	          Dr=SI*condt*(JIp+JI1);
        	          Dz=SI2*condz*JI;
                      %         Dz=SI2*condz*JI;

                      Mdu(1,icont_int)=Dr+Dz;
                    end  %jj
                end   %ii
                ML(icont,:)=Mdu;
                %      'qui ML', keyboard
            end  %jp
        end   %ip
    else %iKcost
        ML=MLm{icnu};
    end

    iML=inv(ML);
    %   'qui inv', keyboard
    if nua==0
        iMLt{icnu}=iML/(2*pi);
    else
        iMLt{icnu}=iML/pi;
    end
    %      'fine calcolo 1 azimut', keyboard
end %azimut
% save tutto1

% ' exist iMLt ?', keyboard

if exist('iMLt')==0
    iMLt=iMLt_sa;
end


fr0=exp(-(ro_out/10).^2);
fi=find(ro_out>20);
fr0(fi)=0;

nVa=2;
if icir==1
    El=0;
else
    El=El0;
end
fr=[];
rt=[];


for kfi=1:NPfi+1
    ri0=Ri0*(1-El*cos(nVa*fian(kfi)));
    fr0=exp(-(ro_out/ri0).^2);
    fi=find(ro_out>20);
    fr0(fi)=0;
    fr=[fr fr0];
    rt=[rt ro_out];
end


fz=exp(-abs(zz-z_act)/2);
fi=find(zz>Lzd | zz<Tbuf);
fz(fi)=0;

[xfm,yfm]=meshgrid(ro_out,zz);




T_VCSEL=sum(Tvet);
z_dd=zf(end);
%z_ddT=zz(end)-z_dd;   % old no metallo
z_ddT=zz(end)-z_dd-Tmetallo;
Tr_buf=zf(end)-T_VCSEL;

Tbuf_shift=zz(end)-z_dd;

%fi_dd=find(zz>=z_ddT);     % old no metallo
%z_ddT=0;
%  'qui terl', keyboard

if abs(zz(end)-Tbuf)>.5
    fi_dd=find(zz>=z_ddT);
else
    fi_dd=1:length(zz);
end
zTdd=zz(fi_dd)-z_ddT;
% zTdd=zz(fi_dd)-Tbuf;
fir=find(xx<=rho(end));
xTdd=xx(fir);
z=z_dd-fliplr(zf);
z=(zf);
%qq=fliplr(Joule+Rec_nr+OptAbs);
%'qui ter', keyboard
%[roM,zM]=meshgrid(rho,z+dZ_ther_dd);
[roM,zM]=meshgrid(rho,z);
[xxd,zzd]=meshgrid(xx,zz);
z1=zzd(:,1);
dzN=[0; diff(z1)];
x1=xxd(1,:);
xdxN=x1.*[0 diff(x1)];
xdxO=xdx;
dzO=dz;
iRAD_spalmato=mode.RAD_spalmato;
if iRAD_spalmato==1
    % 'entro spal', keyboard
    xdx=xdxN;
    dz=dzN;
    %   qtot=Rec_RAD;
    Ra=reshape(mode.alpha,mesh.nny,mesh.nnx);
    RaV=interp2(roM,zM,Ra(2:end-1,:),xTdd,zTdd','linear',0).';
    qtotA=interp2(roM,zM,Rec_RAD',xTdd,zTdd','linear',0).';
    qtot=zeros(length(xx),length(zz));
    qtot(fir,fi_dd)=qtotA;

    %    'qui terl', keyboard
    Rtot=zeros(length(xx),length(zz));
    Rtot(fir,fi_dd)=RaV;
    fiV=find(z1>Z_sub);
    Rint=xdx*Rtot(:,fiV)*dz(fiV);
    RtoT=zeros(size(Rtot));
    RtoT(:,fiV)=Rtot(:,fiV)/Rint;

    qrad1=xdx*qtot(:,fiV)*dz(fiV);
    qradz=qtot(:,fiV)*dz(fiV);
    LZmed=sum(dz(fiV));
    QradMed=qradz/LZmed;
    FZ=zeros(size(dz'));
    FZ(fiV)=1;
    Qver0=xdx*QradMed;
    qtot=QradMed*FZ;
    QtInt=xdx*qtot*dz;

    AlInt=xdx*RtoT*dz;

    Rspalmatod=RtoT.*qtot;
    Qver=xdx*Rspalmatod*dz;
    Rspalmato=Rspalmatod*qrad1/Qver;
    qtot_RAD=Rspalmato;
    %  qq=(Joule+Rec_Cap+Rec_srhAu+OptAbs+Tho);
    qq=(Joule+Rec_Cap+Rec_srhAu+OptAbs);
    %  qq=(Joule+Rec_Cap+Rec_srhAu);
    xdx=xdxO;
    dz=dzO;
else
    qtot_RAD=0;
    %  qq=(Joule+Rec_Cap+Rec_srhAu+OptAbs+Tho+Rec_RAD);
    qq=(Joule+Rec_Cap+Rec_srhAu+OptAbs+Rec_RAD);
    %qq=(Joule+Rec_Cap+Rec_srhAu+Rec_RAD);
    %'passo', keyboard
end
%'passo', keyboard

fiN=find(isnan(qq)==1);
qq(fiN)=0;
% [ZZ,XX]=meshgrid(zTdd,xTdd);
% figure,surf(XX,ZZ,qtotV0)
qtotV0=interp2(roM,zM,(qq).',xTdd,zTdd','linear',0).';
qtot=zeros(length(xx),length(zz));
qtot(fir,fi_dd)=qtotV0;
qtot=qtot+qtot_RAD;

ivecchioBuf=0;
if ivecchioBuf==1
    ' qui controllo termico 0', keyboard
    % Joule heating nel substrato
    fizbuf=find(zz<Tbuf);
    zbuf=zz(fizbuf);
    roBmin=ro_met;
    BuBrod=1.1; % fattore di allargamento corrente nel substrato (bassissimo)
    roBmax=BuBrod*roBmin;
    r_z=roBmax-zbuf*(roBmax-roBmin)/Tbuf;
    Area=pi*roh(1).^2;
    rocond_buf=mode.Zmat/Tbuf*Area;
    Cur=mode.ii_dd(end);
    qbuf=rocond_buf*(Cur./(pi*r_z.^2)).^2;
    qvb=zeros(length(xx),length(fizbuf));
    for k=fizbuf
        rz=r_z(k);
        firo=find(xx<rz);
        qvb(firo,k)=qbuf(k);
    end

    qtot(:,fizbuf)=qvb;
end
%corbuf=
%qBuf=


%   ' qui controllo termico', keyboard

% [zM,roM]=meshgrid(z,ro);
% qtot=interp2(roM,zM,qq,xTdd,zTdd');
% cambio grigliato

% 'qtot', keyboard

if itras==10
    qm=mean(qtot);
    qtotv=qtot;
    qtot=repmat(qm,length(ro_out),1);

    % 'qtot', keyboard
end

Rp=R';
for kd=1:J
    Rdx(kd,:)=Rp(kd,:).*xdx;
end

%    'qui New', keyboard
%dz=diff(ztot);
%dz=[dz dz(end)];
%dz([1 end])=dz([1 end])/2;

Pqs=qtot*(diag(dz)*S.');
P=Rdx*Pqs;


cijL=[];
for inu=1:length(nu_vett)
    Ji=Js(inu);
    Jsh=sum(Js(1:inu-1));
    QL=reshape(P((1:Ji)+Jsh,:)',Ji*L,1);
    cijL=[cijL; iMLt{inu}*QL];
end
CoL=reshape(cijL,L,J)';

ichQ=0;
if ichQ==1
    MMI=sum(MI,3);
    MMIB=sum(MIB,3);
    Q_ric=Ru*inv(MMIB)*P*inv(MMI)*S/2/pi;
    figure, plot(zz,Q_ric,zz,qtot,'.'), pausak
end

T3D1=Ru*CoL*S;

if icir==0
    T3D=reshape(T3D1,length(xx),NPfi+1,length(zz));
    T3Dt=squeeze(mean(T3D,2));
else
    T3D=T3D1;
    T3Dt=T3D1;
end

if iold==0
    Ti=((T_Peltier+T3Dt)/Temp0).^Exp_Temp0;
elseif iold==2
    Et0=-1.3
    Et=Et0+2e-3*(T3Dt+T_Peltier-Temp0);
    %Et=-1.3+4e-3*Tprec;
    Ti=((T_Peltier+T3Dt)/Temp0).^Et;
elseif iold==1
    Ti=(1+T3Dt/Temp0).^Exp_Temp0;
elseif iold==-1
    Ti=1;
end

%if iold==0
% Ti=((T_Peltier+T3Dt)/Temp0).^Exp_Temp0;
%else
% Ti=(1+T3Dt/Temp0).^Exp_Temp0;
%end

%Ti=(1+T3Dt/Temp0).^Exp_Temp0;
Fcondt0=Fcondt;
Fcondz0=Fcondz;

T3D0=T3Dt;
Tprec0=T3Dt;
eRmin=1e-3;

%'iplotTer', keyboard
iplotTerm=0;
if (isfield(mode,'iplotTerm')==1 & IPLOT>0 )
    iplotTerm=mode.iplotTerm;
end
col='mcrgw';
[du,fiat]=min(abs(zz-z_act));

if (iplotTerm & length(Tprec)>1)

    h=figure(38992);clf
    set(h,'pos',[65 66 736 365])
    subplot(1,3,1)
    grid on
    hold on
    box on
    plot(xx,Tprec(:,fiat),'b','linewidth',2)
    ylabel('Temp rise (K)')
    xlabel('\rho (um)')
    a=axis;
    a(2)=2*ro_mesa;
    a(4)=a(4)*1.2;
    axis(a)
    subplot(1,3,2)
    grid on
    hold on
    box on
    plot(zz,Tprec(1,:),'b','linewidth',2)
    a=axis;
    xlabel('z (um)')

    a(1)=Tbuf*.95;
    a(2)=zz(end);
    a(4)=a(4)*1.2;
    axis(a)
    subplot(1,3,3)
    grid on
    hold on
    box on
    plot(zz,Tprec(1,:),'b','linewidth',2)
    xlabel('z (um)')
    % pausak
    drawnow
end  %iplotTerm


if iKcost==0
    % hold on
    % 'qui prima di loop', keyboard
    if icir==0
        T3D=reshape(T3D1,length(xx),NPfi+1,length(zz));
        T3Dt=squeeze(mean(T3D,2));
        Tat=squeeze(T3D(:,:,fiat));
    else
        T3Dt=T3D1;
        Tat=T3Dt(:,fiat);
    end
    %  'qui temnp', keyboard

    PTherm=2*pi*xdxN*qtot*dzN*1000;
    % 'prima di termici', keyboard

    if mode.ioldTemp==-1
        fattore_correttivo=1
    else
        % fattore_correttivo=mode.PDiss(end)/PTherm;
        fattore_correttivo=mode.PDissPred(end)/PTherm
    end
    %Tdd=Tdd*fattore_correttivo;
    if max(max(Tprec))>7
        %'qui tempo', keyboard
    end
    for kit=1:10

        Fcondt=Fcondt0.*Ti';
        Fcondz=Fcondz0.*Ti';

        subTEMP

        if icir==0
            T3D=reshape(T3D1,length(xx),NPfi+1,length(zz));
            T3Dt=squeeze(mean(T3D,2));
            Tat=squeeze(T3D(:,:,fiat));
        else
            T3Dt=T3D1;
            Tat=T3Dt(:,fiat);
        end

        if iold==0
            Ti=((T_Peltier+T3Dt)/Temp0).^Exp_Temp0;
        elseif iold==2
            Et0=-1.3
            Et=Et0+2e-3*(T3Dt+T_Peltier-Temp0);
            %Et=-1.3+4e-3*Tprec;
            Ti=((T_Peltier+T3Dt)/Temp0).^Et;
        elseif iold==1
            Ti=(1+T3Dt/Temp0).^Exp_Temp0;
        elseif iold==-1
            Ti=1;
        end


        %if iold==0
        % Ti=((T_Peltier+T3Dt)/Temp0).^Exp_Temp0;
        %else
        % Ti=(1+T3Dt/Temp0).^Exp_Temp0;
        %end

        %Ti=((T_Peltier+T3Dt)/Temp0).^Exp_Temp0;
        %   Ti=(1+T3Dt/Temp0).^Exp_Temp0;
        %   su=sum(sum(Tprec));
        %   sui=sum(sum(T3Dt));
        su=max(max(Tprec));
        sui=max(max(T3Dt));
        if su>10
            er0=abs(1-su/sui)
        else
            er0=0;
        end
        er(kit)=er0;
        ishow=0;
        %   'qui iter new', keyboard

        Tprec=T3Dt;
        if ishow==1
            'iter =', kit
            'errore relativo Temp ',

            figure(h);
            subplot(121)
            plot(xx,Tprec(:,fiat),col(kit)),
            a=axis;
            a(2)=2*ro_mesa;
            axis(a)
            subplot(122)
            plot(zz,Tprec(1,:),col(kit)),
            a=axis;
            a(1)=Tbuf*.95;
            a(2)=zz(end);
            axis(a)
            pausak
        end

        %keyboard
        if er0<=eRmin
            break
        end

    end
else
    if icir==0
        T3D=reshape(T3D1,length(xx),NPfi+1,length(zz));
        T3Dt=squeeze(mean(T3D,2));
        Tat=squeeze(T3D(:,:,fiat));
    else
        T3Dt=T3D1;
        Tat=T3Dt(:,fiat);
    end

end  %iKcost

if ifig==1 & length(Tprec)>1

    [du,fiat]=min(abs(zz-z_act));
    figure,
    subplot(211)
    plot(xx,Tprec0(:,fiat),xx,Tprec(:,fiat)),
    a=axis;
    a(2)=2*ro_mesa;
    axis(a)
    subplot(212)
    plot(zz,Tprec0(1,:),zz,Tprec(1,:)),
    a=axis;
    a(1)=Tbuf*.95;
    a(2)=zz(end);
    axis(a)
    pausak
end

toc
if ifig==1
    figure
    %      surf(xfm',yfm',squeeze(T3D(:,1,:))), shading interp, view(2), colorbar,
    surf(xfm',yfm',T3Dt), shading interp, view(2), colorbar,
    %      title([' Cor = ',num2str(Cor)])
    pausak




    figure, plot(ro_out,Tat,'.-')
    title(' Temp zona attiva ')
    pausak

end


Tbuf_aggiunto=max(zz)-max(z);

rho(1)=min(min(xxd));

%'Termico', keyboard

Tdd=interp2(xxd,zzd,T3Dt',rho,Tbuf_aggiunto+z');
%  Tdd=flipud(Tdd);
Tdd=[Tdd(1,:);Tdd;Tdd(end,:)]; % to restore the contact



Contributi_Termici     % T_Contributi is computed

for iTc=1:length(T_Contributi)
    T_Contributi{iTc}=interp2(xxd,zzd,T_Contributi{iTc}',rho,Tbuf_aggiunto+z');
    %  Tdd=flipud(Tdd);
    T_Contributi{iTc}=[T_Contributi{iTc}(1,:); T_Contributi{iTc}; T_Contributi{iTc}(end,:)]; % to restore the contact
end


K = interp2(xxd,zzd,Fcondz,rho,Tbuf_aggiunto+z');
K = [K(1,:);K;K(end,:)];
K = reshape(K,1,mesh.nn);

TotalHeat = interp2(xxd,zzd,qtot',rho,Tbuf_aggiunto+z');
TotalHeat = [TotalHeat(1,:);TotalHeat;TotalHeat(end,:)];
TotalHeat = reshape(TotalHeat,1,mesh.nn);

if iploCond==1
    figure, surf(xxc,zzc,log10(Fcondt)), view(2), shading flat,
    title(' Transverse Thermal conducibility')
    a=axis;
    a(3)=Tbuf-2;
    a(4)=max(max(zzc))+1;
    axis(a)
    colorbar
    caxis([-7 -6]+2)
    pausak
end
%pausak
%'fine termico', keyboard
