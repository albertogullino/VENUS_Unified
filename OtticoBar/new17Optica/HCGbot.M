function hcgbot(ini,dret,drel,dox,dsub,period,DC,fi_add,pun,ipro,iprec,ifp,...
 alim,iLP,ianyv,dk0,dec_in,dpm_in,H,r_lens,Ndisc,Lbuf,dmet,Ps,na0,amod,Pf,pasnu,fun_wolu)
colordef('black')
%clear global
close all
%path(path,'D:\USERS\PIERLUIGI\MVCSEL\new03')
%path(path,'D:\USERS\PIERLUIGI\MVCSEL\new02')
%ifp=-4;
%ifp=-10;
%ipro=1;
isemp=1;
chr='or';
iint=1;   % cerca intervallo
if ifp~=-4
 iint=0;
end
iint=0;
ipolar=1;
iopt=0;

%ini
%keyboard
if isfield(Ps,'Nmiru')==0
 Nmiru=20;
else
 Nmiru=      Ps.Nmiru;
end
if isfield(Ps,'dat')==0
 dat=dox;
else
 dat=      Ps.dat;
end
if isfield(Ps,'Delret')==0
 Delret=0;
else
 Delret=      Ps.Delret;
end
%'ikr', keyboard
if isfield(Ps,'ikr')==0
 ikr=0;
else
 ikr=      Ps.ikr;
end
if isfield(Ps,'tolk')==0
 tolk(1)=.05;
 tolk(2)=.0;
else
 tolk=      Ps.tolk;
end
if isfield(Ps,'farm')==0
 farm=1;
else
 farm=      Ps.farm;
end
t0=clock;

renom=[1 4 3 2 5 6];
varvel=dret;
varlen=drel;


%Ev_or_Od='Both';
%Ev_or_Od='Odd';
%Ev_or_Od='Even';
%Ev_or_Od='[5]';
%'amod', keyboard
if amod>0
 Ev_or_Od=num2str(amod);
else
 if abs(amod)==1
  Ev_or_Od='Even';
 elseif abs(amod)==2
  Ev_or_Od='Odd';
 elseif abs(amod)==3
  Ev_or_Od='Both';
 end


end


 shf0=0.2;
 FRint0=-.1;
 FRint=.1;
 NPfreq=3;


 fin0='';
 np_miru=-100;
 np_mir=Nmiru;

 dpm=37.4;
 dpm1=0;
 fil_str=[ini,'.str'];


if nargin<14
 iLP=0;
end

%nK_dis=30;
%nK_dis=20;

nK_dis=Ps.Nk;

if nargin<13
 alim=.25;
end
%dk=.004;


if iprec==5
 nuac_f1=9;
 nuac_f0=7;
% dk=.005;
 dk=.006;
elseif iprec==4
 nuac_f1=8;
 nuac_f0=5;
 dk=.006;
elseif iprec==3
 nuac_f1=6;
 nuac_f0=4;
 dk=.007;
elseif iprec==2
 nuac_f1=5;
 nuac_f0=3;
 dk=.007;
 dk=.01;
elseif iprec==1
 nuac_f1=4;
 nuac_f0=2;
 dk=.008;
% dk=.012;
elseif iprec==0
 nuac_f1=4;
 nuac_f0=2;
 dk=.01;
end
%'vedo', keyboard
% nuac_f1=6;
tich=num2str(pun);
if pun<0
 pun=0;
 tich='00';
end

if pun==20
if iprec==5
 dk=.002;
elseif iprec==4
 dk=.003;
elseif iprec==3
 dk=.004;
elseif iprec==2
 dk=.005;
elseif iprec==1
 dk=.006;
end
end

nucompi0=na0;

%'na0', keyboard



if nargin<15
 iany=1;
end
if nargin>=16
 if dk0>0
  dk=dk0;
 end
end
shfom =0.5;
shfom =0.;

unpun=0;


%varlen=[0];
%varvel=[50];


%DC=0.5;
%pun=[0 1 2];
%pun=[0];

tistru=2;    %0 rel, 1 gratrel, 2 gratrel all (reticolo ovunque), 3 rel+gratrel, 4 rel+grat_ext, 5 no-rilievo(solo reticolo), 6 come 4 ma senza reticolo dentro


%for isiz=[1 2]
%for isiz=[1 ]

for tistrup=pun
%tistrup, keyboard
 if tistrup==0
  tistru=tistrup;
%  ' tistru ', keyboard
  if alim==0
   alim=0.28;
  end
%  ikr=0;
  nucompi=0;
  shfom =shf0;
 else
%  ikr=2;
  nucompi=nucompi0;
  shfom =0.5;
  shfom =shf0;
  tistru=renom(tistrup);
 end

 %DC=0.5;
% if tistru==4
%  DC=0.2;
% end
% fi_sa=[fin0,num2str(tistrup),num2str(isiz)]
% fi_sa=[fin0,num2str(tistrup),chr(isiz)]
% fi_sa=[fin0,chr(isiz),num2str(tistrup)]
% fi_sa=[fin0,fi_add,tich]
 fi_sa=fi_add;
% keyboard
 set_mes


if exist('dpm_in')==1
 if exist('dpm')==0
  dpm=dpm_in;
 end
end

 if nargin>=17
  if dec_in>0
   dpm=dpm-dec_in;
  end
 end


 if ipro==1
  'check loops ', keyboard
 end

if exist('Pf')
 if length(Pf)>0
 shfom=Pf.sh;
 FRint0=Pf.fi;
 FRint=Pf.fu;
 NPfreq=Pf.N;
 ipolar=Pf.ipolar;
 end
end

P.pasnu=pasnu;
P.ipolar=ipolar;
% save data iLP fi_sa ifp alim nK_dis numodiacc ivestru tolk dk ipro Dlam_mod Ev_or_Od Dreti Dreli nucompi ikr iany dsub tistru DC period isem
 P.fin0=fin0;
% ' qui ',  keyboard

 P.farm=farm;
 P.FRint=FRint;
 P.FRint0=FRint0;
 P.iLP=iLP;
 P.shfom =shfom;
 P.fi_sa=fi_sa;
 P.ifp=      ifp;
 P.alim=     alim;
 P.nK_dis=   nK_dis;
 P.numodiacc=numodiacc;
 P.tolk=     tolk;
 P.dk=       dk;
 P.ipro=     ipro;
 P.Pf=Pf;
% P.Dlam_mod= Dlam_mod;
 P.NPfreq=   NPfreq;
 P.Ev_or_Od= Ev_or_Od;
 P.Dreti=    Dreti;
 P.Dreli=    Dreli;
 P.nucompi=  nucompi;
 P.ikr=      ikr;
  P.iany=     ianyv(1);
 if length(ianyv)==2
  P.ianys=     ianyv(2);
 else
  P.ianys=     0;
 end
 P.dsub=     dsub;
 P.tistru=   tistru;
 P.DC=       DC;
 P.Del_su=       Delret;
 P.period=   period;
 P.dpm=      dpm;
 P.dpm1=     dpm1;
 P.fil_str=  fil_str;
 P.np_mir=   np_mir;
 P.np_miru=   np_miru;
 P.H=   H;
 P.r=   r_lens;
 P.Ndisc=   Ndisc;
 P.Lbuf=   Lbuf;
 P.dmet=   dmet;
 P.dox=   dox;
 P.att=   dat;
 %'Ps',   keyboard
 P.Ps=Ps;


% gen_wol(P,iint)

% gen_woluBu(P,iint)
% 'doxe', keyboard

 global fun_wolu

if length(fun_wolu)==0
 fun_wolu='gen_woluB';
end
 
 %'doxe', keyboard

   fu=[fun_wolu,'(P,iint)'];
   eval(fu)



end
%end % isiz

   tim=etime(clock,t0);
   ore=fix(tim/3600);
   minuti=fix((tim-ore*3600)/60);
   secondi=(tim-ore*3600-minuti*60);
   format short
   STRI=[' TOTAL ',ini,' STRUCTURE SIMULATION time (hours, minutes, seconds) = '];
   disp(STRI);
   disp([ ore minuti secondi]);
