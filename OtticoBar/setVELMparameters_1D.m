
%-- macro per settare i parametri di VELM

% colordef black

calop=struct(...
'Dlam',[],...
'mod_max',[],...
'ro_in',[],...
'zeta',[],...
'epsT',[],...
'N',[],...
'Nref',[],...
'numodiacc',[],...
'alim',[],...
'dk1',[],...
'par_in',[],...
'Ev_or_Od',[],...
'ipolar',[],...
'pasnu',[],...
'irid_bas',[],...
'ispeed',[],...
'iPERD_vol',[],...
'iplan',[],...
'del_n_ag',[],...
'pimra',[],...
'ianti_gui',[],...
'box',[]...
);

% parametri da non modificare
ispeed=0;
irid_bas=1;
iPERD_vol=1;
iplan=0;
calop.ro_in=0;
calop.zeta=0;
calop.epsT=0;
calop.N=0;

% calop.NPlam=5;
% calop.dlam=.01;
calop.NPlam=VelmOptions.NPlam;
calop.dlam= VelmOptions.dlam;


%'setopt_dd', keyboard

calop.mod_max  =mod_max;
calop.numodiacc=numodiacc;
calop.alim     =alim;
calop.dk1      =dk1;
calop.Ev_or_Od =Ev_or_Od;
calop.ipolar   =ipolar;
calop.pasnu    =pasnu;
calop.ispeed   =ispeed;
calop.irid_bas   =irid_bas;
calop.Dlam     =Dlam;
calop.iPERD_vol     =iPERD_vol;
calop.iplan=iplan;
%ianti_gui
%'setopt_dd', keyboard


 ifpO=-4;   % set to -4 to overwrite the .par file derived from the str
 Npar=gest_opt(fil_str,ifpO);
 
 fip=findstr(fil_str,'.');
 nomestr=fil_str(1:fip-1);
% ' dopo gest_opt', keyboard

 [cPAR,cPver]=gest_par(nomestr,ifpO);
 Npara=length(cPver);
% ' dopo gest_par', keyboard

%par_in=cell(Npar,1);

pin=[];
REM=cPAR;
PARV=par_in;
if Npara~=length(PARV)
 ' parametri non corretti !   ',
 ' far coincidere PARV (in parv.m) coi parametri @ in ordine   ', keyboard
 return
end

for ip=1:Npara
 pin{ip,1}=PARV{ip};
 [ch_p{ip,1},REM]=strtok(REM);
end


%if length(pin)>0
% par_in=loop_ass(par_in,cPver,cPAR,pin,ch_p,0);
 calop.par_in=par_in;
%end

%' qui set_opt', keyboard


global pimra ipim ifpdis
global Nr_dyn c_light tauN0 hbar

hplanck=6.626070040e-34; % Plank constant, J*s
hbar=hplanck/(2.*pi);
c_light=2.99792458e+8; % light velocity, m/s
mi=4*pi*1.0e-7; % H/m
eps0=1/(mi*c_light^2); % F/cm

tauN0=1;

Nr_dyn=NPx1;    % discret. radiale corrente e temperature

ifpdis=ifp_Opt;
ipim=0;
%pimra=-0.01;
ipim=1;
pimra=0;
% c0=3e-1; %in micron/ns
% ra=real(rqw)+j*pimra*NQW;
% cconv=k0*c0/rr;
% g0=pimra*cconv*NQW;

c0=c_light*100;
k0=2*pi/lambda*1e4;
cconv_velm=2*k0;
%'cconv', keyboard

MODC0=0;

Ps.iord_long=1;

iBEL=100 ;   % reticoli
%'qui velma', keyboard
iexis=isfield(mode1,'TmVelm');
if (isfield(mesh,'DeltaTvelm') & iexis==1 & isfield(mode,'Gmod'))
 dT_nuovo=max(max(mesh.DeltaTvelm));
 Tv=mode1.TmVelm;
 Lv=mode1.LamVelm;
%'qui velma', keyboard
%mode.mintempVELM
 %dTmin0=mode.mintempVELM;
 %dTmin=5;
 %fiT=find(Tv>dTmin0);
 gr=1;
 if isfield(mode,'DT0')
  DT0=mode.DT0;
 else 
  DT0=10;
 end 
% Tv
% mode.Gmod

 maGmod=max(mode.Gmod,[],1);
   if Ps.Pf.ipolar==3  
%    'qui setVelm', keyboard
   end  
% 'qui set', keyboard   
 if maGmod(end-2)>0 & length(Lv)>gr+1 &  Tv(end)>DT0
 %if Tv(end)>DT0 & mode.Gmod(end-2)>0
 %if length(Tv)>5 & mode.Gmod(end-2)>0
   if Ps.Pf.ipolar==3  
%    'qui set', keyboard
   end 
  fiTp=sort(length(Lv)-[0:gr+1]);
  LvT=Lv(fiTp);
  TvT=Tv(fiTp);

  coLT=polyfit(TvT,LvT,gr);
  Lvf=polyval(coLT,TvT);
   Lvf_nuovo=polyval(coLT,dT_nuovo);
%       figure(63), plot(Tv,Lv,TvT,Lvf,'.',dT_nuovo,Lvf_nuovo,'ro'), pausak
%' qui set', keyboard
   if mode.verbVELM>=1
    figure(63), plot(Tv,Lv,TvT,Lvf,'.',dT_nuovo,Lvf_nuovo,'ro'), pausak
%    figure, plot(Tv,Lv-Lvf), pausak
   else
%      figure(63), plot(Tv,Lv,TvT,Lvf,'.',dT_nuovo,Lvf_nuovo,'ro'), 
   end
   TvT=[TvT dT_nuovo];
   Lvf=[Lvf Lvf_nuovo];
   dL=diff(Lvf);
   ddL=dL(1:end-1).*dL(2:end);
%   if prod(ddL)>0
    iBEL=99;

   Ps.lambda_Temper=Lvf_nuovo;
%   Dlam(1:2)=VelmOptions.Dlam(1:2)-mean(VelmOptions.Dlam(1:2))/3;
   
   if length(mode.vlambda)>1
    Dlam(1)=-abs(Dlam(5));
    DDL=(max(mode.vlambda)-min(mode.vlambda))    ;
    Dlam(2)=1.3*DDL;
   else
    Dlam(1)=-abs(Dlam(5));
    Dlam(2)=abs(Dlam(5));
   end
   if mode.nmodes==1
    Dlam(3)=3;
    if length(Dlam)>=6
     Dlam(3)=Dlam(6);
    end
   else
    Dlam(3)=max([ 3 ceil(diff(Dlam(1:2))/abs(2*Dlam(5)))]);
   end 
   calop.Dlam     =Dlam;   
   end  
   %'qui', keyboard

% end  
end
%Dlam(1:2)=VelmOptions.Dlam(1:2);
%iBEL=101;
%   calop.Dlam     =Dlam; 
%   'qui', keyboard

Ps.iBEL=iBEL;



Ps.dissfun='diss_new';

if isfield(VelmOptions,'dissfun')
 Ps.dissfun=VelmOptions.dissfun;
end

if(verbVELM<=0)
    % Tibaldi: questa riga permette di tagliare molti messaggi
    Ps.ifpstop=0; 
else
    Ps.ifpstop=1;
end
%'qui Az', keyboard
if isfield(mode1,'campoz')
if mode1.campoz==1
'qui Az', keyboard

% settings per plot campo ovunque
Ps.isav_Az=1;
Ps.izoom=1;
Ps.Lizi=20/1000;    % step longitudinale in micron   %
end
end
calop.Ps=Ps;
%'FERP', keyboard