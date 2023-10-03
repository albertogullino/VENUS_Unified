  global isa_file nfile_sa  lpvc
  global lp1 lp2
  if length(isa_file)==0
   isa_file=0;
  end

  if isa_file==1
    nomsave=[nfile_sa,num2str(iLP),'_',num2str(mm)];
  elseif isa_file==2
    nomsave=[nfile_sa,lpvc,'_',num2str(iLP),'_',num2str(mm)];
  end
%           nomsave
%  ' sub_old', keyboard


%  isaMo=1;


  T0=Temp;
  KKsav=0;
  iraffi=0;
  Pol.Ex=[];
  Dla_out=[];

  clear Anuvet Acavet Amevet Gvet alvet Fint Kvet
  clear Anuvet1 Acavet1 Amevet1 Gvet1 alvet1
  clear Gvet2 alvet2
%  clear EDout EDc Gsov Fsov
%  clear Nazim Nrad tyE Wv czv fqwv fa4v Lefv Tefv Tefiv
  clear Azvet Azvet1 Azvetf Azvetf1 Azvetf2 Azvet2

global ispeed sl_dev
global gou1 aou1 gou2 aou2 fou1

%'ispeed', keyboard
if length(ispeed)==0
 ispeed=1;
end




  nube=mm;

if mm~=0 &  mm~=1 & numodiacc>0
 ' subold non adatta ', keyboard
end
  nubesi=nube-pasnu*numodiacc;
  nubesi=nube;
%  nubesi=nube-pasnu*numodiacc;
%  if nubesi<-1
%     if mm==1
%      nubesi=-1;
%    end
%    if mm==0
%     nubesi=-2;
%   end
% end
%  disp(' provvisorio nubewsi'), keyboard
%  nubesi=2;
  nubesu=nube+pasnu*numodiacc;
 if pasnu==1
  nubesu=nubesi+numodiacc;
 end
 
 if ipolar~=0
  numodi=numodiacc+1;
 else
  numodi=(numodiacc+1)*2;
%  numodi=((nubesu-nubesi)/pasnu+1);
end

numodi
%'numodi', keyboard
%  fre_camp=linspace(Frisi,Frisu,Ndisp);
%  fre_camp=linspace(Frisu,Frisi,Ndisp);
%'fre_camp', keyboard

  eint=alim;


  if iLP==1
   nni=numodi;
  else
   nni=2*numodi;
  end

  pes0=1;

%  disp('sub'), keyboard
imem=1;   %=1 velocizza le matrici diagonali
%imem=0;   %=1 velocizza le matrici diagonali
%'ipolar', keyboard

if ipolar==2
  pvet=[1 -1];
else
  pvet=ipolar;
end
pvet00=pvet;

%load lam

for ifr=1:length(fre_camp)

 if ifr==2
%  ' controllo ispeed', keyboard
 end

 
  freq=fre_camp(ifr);

%  if ifast==0
     kcav=kcav0*(1+freq);
%  end   
  
  if ivfre==1
   n_update
  end
%  if ifp~=-4
  [ifr freq length(fre_camp)]
%  end

  %kcav=kcav0*(1+freq);

%      disp(' prima acc0 '), keyboard

 ' load stop '
  load stop
  if istop==1
   keyboard
  end

   tempo0=clock;
%'acc0', keyboard
   acc0
%   if ifp~=-4
%   unafreq=etime(clock,tempo0)
   if exist('ifp_post')
    'dopo acc0 ', keyboard
   end
%   end
if ifr==2
% ifp=-12
end 
% N1i=-inv(N1);
% Mtot=N1i*N2;
% mapab(Mtot)

%      disp(' dopo acc0  in sub_mix')
%      keyboard


   Fint(ifr)=freq;

%   if ipolar==-1 | ipolar==2
%   if ipolar==1 | ipolar==2
   if length(pvet)>=1

    sAz=size(Anz1);
    if length(sAz)==3
     Azvet1(1:sAz(1),1:sAz(2),1:sAz(3),ifr)=Anz1;
     Azvetf1(1:sAz(1),1:sAz(2),1:sAz(3),ifr)=Anzf1;
    else
     Azvet1(1:sAz(1),1:sAz(2),ifr)=Anz1;
     Azvetf1(1:sAz(1),1:sAz(2),ifr)=Anzf1;
    end
    Gvet1(1:sAz(2),ifr)=Gvav1;
    sk=size(KK);
    Kvet1(1:max(sk),ifr)=KK;
    alvet1(1:sAz(2),ifr)=alphavv1;

   end

%   if ipolar==1 | ipolar==2 | ipolar==0
%   if ipolar==-1 | ipolar==2
   if length(pvet)==2

%   ' qui size ', keyboard
    sAz=size(Anz2);
    if length(sAz)==3
     Azvet2(1:sAz(1),1:sAz(2),1:sAz(3),ifr)=Anz2;
     Azvetf2(1:sAz(1),1:sAz(2),1:sAz(3),ifr)=Anzf2;
    else
     Azvet2(1:sAz(1),1:sAz(2),ifr)=Anz2;
     Azvetf2(1:sAz(1),1:sAz(2),ifr)=Anzf2;
    end
    Gvet2(1:sAz(2),ifr)=Gvav2;
    alvet2(1:sAz(2),ifr)=alphavv2;

   end


end % while ifr
%' sub key', keyboard

%if ifast==1
% clear Tstorbm  Tstortm
% clear Tstorbp  Tstortp
%end

  global isa_file nfile_sa  lpvc
  global lp1 lp2
%  ' glob', keyboard
  if length(isa_file)==0
   isa_file=0;
  end
  if isa_file==0
%    save sacont
  elseif isa_file==1
    nomsave=[nfile_sa(1:end-1),num2str(mm),'_',num2str(lp1),num2str(lp2)];
    eval(['save ',nomsave])
  elseif isa_file==2
    nomsave=[nfile_sa,lpvc,'_',num2str(iLP),'_',num2str(mm)];
    eval(['save ',nomsave])
  end

global istop_fie
if istop_fie==1
 ifpsalvo=ifp;
 ifp=-10
end 


  ifps=ifp; if ifp>=-2; ifp=-1; end

%ifp=-10
%   'subold', keyboard
%   'subold', keyboard
%global  uL
   Ppol.fatqw=fatqw;
%   Ppol.gpla=gpla;
   if exist('uL')
    Ppol.uL=uL;
   end 

  if ipolar==2
   pola=pvet(1);
   icpo=1;
   Azvet=Azvet1;
   Azvetf=Azvetf1;
   Gvet=Gvet1;
   alvet=alvet1;
   Ppol.Az1=Azvet1;
   Ppol.Az2=Azvet2;
   Ppol.Azf1=Azvetf1;
   Ppol.Azf2=Azvetf2;
   Ppol.G1=Gvet1;
   Ppol.G2=Gvet2;
   Ppol.A1=alvet1;
   Ppol.A2=alvet2;

  else
    pola=ipolar;
    icpo=1;
    %'passo sub', keyboard 
    Azvet=Azvet1;
    Azvetf=Azvetf1;
    Gvet=Gvet1;
    alvet=alvet1;
    Ppol.Az1=Azvet1;
    Ppol.Azf1=Azvetf1;
    Ppol.G1=Gvet1;
    Ppol.A1=alvet1;
%   pola=ipolar;
%   if ipolar==1 | ipolar==-1
%    Azvet=Azvet1;
%    Azvetf=Azvetf1;
%    Gvet=Gvet1;
%    alvet=alvet1;
%    Ppol.Az1=Azvet1;
%    Ppol.Azf1=Azvetf1;
%    Ppol.G1=Gvet1;
%    Ppol.A1=alvet1;
%   elseif ipolar==-1
%    Azvet=Azvet2;
%    Azvetf=Azvetf2;
%    Gvet=Gvet2;
%    alvet=alvet2;
%    Ppol.Az2=Azvet2;
%    Ppol.Azf2=Azvetf2;
%    Ppol.G2=Gvet2;
%    Ppol.A2=alvet2;
%   end
  end

  if iplan==1 & KK(1)>.05
  ifp=-4;
  end
  isalva=1;
%  ' primo camv_new: ipolar ', keyboard

 if isfield(Ps,'ifpstop')==1
  if Ps.ifpstop==1
   ifpsalvo=ifp;   
   ifp=-10;
  end
 end  
  camv_new


  aou1=aou;
  gou1=gou;
  fou1=fou;



%gsod=gsov;
%fi=find(gsod<0);
%gsod(fi)=1e30;
%[du,imi]=min(gsod);
%lambda_ver=lambda*(1-fsov(imi));

%'qui ver', keyboard

  Ppol.Fint=Fint;
  Ppol.nube=nube;
  Ppol.na0=nucomp;
  Ppol.a0ref=a0ref;
  Ppol.aiat=ar.cor;
  Ppol.istrumix=istrumix;
  Ppol.xroJ=xroJ;
  Ppol.aitot=aitot;
  Ppol.icampi=icampi;
  Ppol.pvet=pvet;
  Ppol.ldap=ldap;
  Ppol.Pus=Pus;
  Ppol.isi=isi;
  Ppol.iFF=iFF;
  Ppol.z0c=z0c;
  Ppol.rr=rr;
  Ppol.nk1max=nk1max;
  Ppol.iztm=iztm;
  Ppol.bvero=bvero;
  Ppol.aou1=aou1;
  Ppol.gou1=gou1;
  Ppol.fou1=fou1;
   Ppol.KK=KK;
   Ppol.KKt=KKt;
   Ppol.KKv=KKv;
   Ppol.numodi=numodi;
   Ppol.pasnu=pasnu;
   Ppol.lbv=lbv;
   Ppol.kcav0=kcav0;
   Ppol.x=x;
   Ppol.fian=fian;
   Ppol.Nx=Nx;
   Ppol.mbv=mbv;
   Ppol.nube=nube;
   Ppol.axtot=axtot;
   Ppol.aytot=aytot;
   Ppol.pola0=pvet;
   Ppol.iLP=iLP;
   Ppol.fimaxi=fimaxi;
   Ppol.lfi_inp=lfi_inp;
   Ppol.ldap=ldap;
   Ppol.Pus=Pus;
   Ppol.isi=isi;
   Ppol.rr=rr;
   Ppol.mbv1=mbv1;
   Ppol.pimu=pimu;
   Ppol.rfu=rfu;
   Ppol.iFFsect=iFFsect;
   Ppol.iFF=iFF;
   Ppol.iFFte=iFFte;
   Ppol.iant=iant;
   Ppol.pes=pes;
   Ppol.z0c=z0c;
   Ppol.lambda=lambda;
   if exist('aax')
   Ppol.aax=aax;
   end
   if exist('Cug')
   Ppol.Cug=Cug;
   end
   Ppol.Litot=Litot;
   Ppol.nitot=nitot;
   Ppol.aitot=aitot;
   Ppol.fmlstot=fmlstot;
   Ppol.par_in=par_in;




%  if (ifp>=-3 | ifp==-10 | ifp==-4) & icampi>0
%  if (ifp>=-3 | ifp==-10 ) & icampi>0 & length(fso)>0
  if exist('aou')
   saou=length(aou);
  end
  if iplan==1 & KK(1)>.05
  ifp=-10;
  end
  if (ifp>=-3 | ifp==-10 ) & saou>0
   vg=vgconv;
   vg1=vg/2;
   fsaf=figure,
   puf=1:size(aou,1);
   lco=lambda*1000;
   zeromo=ones(size(Fsov))*alpha_th;
%c   subplot(211), plot(lco*fou(:,puf),aou(:,puf),lco*fso,zeromo,'wo'), grid;
   subplot(211), plot(lco*fou(puf,:)',aou(puf,:)',lco*Fsov,zeromo,'wo'), grid;
   title(' polar=1 ')
   dista=(max(gou(puf,:))-min(gou(puf,:)))/min(gou(puf,:));
   dista=20;
   if dista>10
    subplot(212), semilogy(lco*fou(puf,:)',gou(puf,:)'/vg1,lco*Fsov,Gsosa/vg,'wo');
   else
    subplot(212), plot(lco*fou(puf,:)',gou(puf,:)'/vg1,lco*Fsov,Gsosa/vg,'wo');
   end
   grid
pausak
   raf=[nomeFs(1:end-4),num2str(iLP),num2str(mm),'P1'];
%   eval(['hgsave(',num2str(fsaf),',''',raf,''');']);
  end

%  if iplan==1 & KK(1)>.05
%   soldi
%   keyboard
%  end

  iscrivo=0;

%  if length(Pol.Ex)>1 & icampi>0
 if icampi>0
  if length(Pol.Ex)>1
   iscrivo=1;
   Ppol.Ex1=Pol.Ex;
   Ppol.Ey1=Pol.Ey;
   Ppol.Cug=Cug;
   Ppol.aax=aax;
   Ppol.XP=XP;
   Ppol.YP=YP;
   if length(find(shailoop==5))>0
    Ppol.cce=cce;
   end
   Ppol.Tfas=Tfas;
   Ppol.Ffas=Ffas;
  end
 end

 if ipolar==2
  icpo=2;
  Azvet=Azvet2;
  Azvetf=Azvetf2;
  Gvet=Gvet2;
  alvet=alvet2;
  pola=pvet(2);
  isalva=1;
%  ' secondo camv_new: ipolar ', keyboard
  camv_new

  aou2=aou;
  gou2=gou;
  fou2=fou;
  Ppol.aou2=aou2;
  Ppol.gou2=gou2;
  Ppol.fou2=fou2;

 if icampi>0
  if length(Pol.Ex)>1 & icampi>0
%  if length(Pol.Ex)>1 & icampi>0
   Ppol.Ex2=Pol.Ex;
   Ppol.Ey2=Pol.Ey;
   if iscrivo==0
    Ppol.Cug=Cug;
    Ppol.aax=aax;
    Ppol.XP=XP;
    Ppol.YP=YP;
    if length(find(shailoop==5))>0
     Ppol.cce=cce;
    end
    Ppol.KK=KK;
    Ppol.KKt=KKt;
    Ppol.numodi=numodi;
    Ppol.pasnu=pasnu;
    Ppol.lbv=lbv;
    Ppol.kcav0=kcav0;
    Ppol.x=x;
    Ppol.fian=fian;
    Ppol.Nx=Nx;
    Ppol.mbv=mbv;
    Ppol.Tfas=Tfas;
    Ppol.Ffas=Ffas;
   end
  end
 end

%  if (ifp>=-3 | ifp==-10 | ifp==-4) & icampi>0
%  if (ifp>=-3 | ifp==-10 ) & icampi>0 & length(fso)>0
  if exist('aou')
   saou=length(aou);
  end
  if (ifp>=-3 | ifp==-20) & saou>0
   figure(fsaf),
   vg=0.5;
   puf=1:length(fso);
   zeromo=ones(size(fso))*alpha_th;
   lco=lambda*1000;
%   hold on
   subplot(211), hold on
% plot(fos,aos,zos,zeros,'wo'), grid;
   plot(fos,aos,zos,zers,'wo',...
   fou(puf,:)'*lco,aou(puf,:)','--',fso*lco,zeromo,'wo'), 
   dista=(max(gou(puf,:))-min(gou(puf,:)))/min(gou(puf,:));
   title([num2str(tyE),'---->   1:X, -1:Y  '])
 %  hold on
%   gos=gou(puf,:)'/vg;
%   zosg=gso/vg;
%   hold on, plot(fos,gos,zos,zosg,'wo',...
%   fou(puf,:)'*lco,gou(puf,:)'/vg,'--',fso*lco,gso/vg,'wo'), grid;
%     plot(fso*lco,gso/vg,'ws',fou(puf,:)'*lco,gou(puf,:)'/vg,'--');
    subplot(212), hold on, semilogy(lco*fou(puf,:)',gou(puf,:)'/vg,'--',lco*fso,gso/vg,'ws'); grid
   ax=axis;
   ax(3)=min(gso/vg)*2/4;
   ax(4)=max(gso/vg)*4/2;
%   axis(ax)
   raf=[nomeFs(1:end-4),num2str(iLP),num2str(mm),'P2'];
   
%   eval(['hgsave(',num2str(fsaf),',''',raf,''');']);
  end

 end
 ifp=ifps;
% 'ispeed sub'
% keyboard

% 'ispeed sub dopo'
% keyboard

%if ifp==-4
% keyboard
% dispfig
% keyboard
%end
%  end
if ifp==-10
'subold'
 keyboard
 keyboard
end

% keyboard
 if icampi>2
  if ifp>-4
   keyboard
  end
 end

  PD.a=zeros(size(Gsov));
  PD.b=zeros(size(Gsov));
  PD.c=zeros(size(Gsov));
  PD.d=zeros(size(Gsov));
  PD.e=zeros(size(Gsov));
  PD.f=zeros(size(Gsov));
  PD.g=zeros(size(Gsov));
  PD.r=zeros(size(Gsov));
  PD.s=zeros(size(Gsov));
  PD.t=zeros(size(Gsov));
  PD.e1=zeros(size(Gsov));
  PD.e2=zeros(size(Gsov));
  PD.f1=zeros(size(Gsov));

%' PD: pausa ', keyboard
 if idyn==1
  PD.a=taut;
  PD.b=czv;
  PD.c=fqwv;
  PD.d=Lefv;
  PD.eu=teff;
  PD.e=tauu;
  PD.f=taub;
  PD.e1=pvol;
  PD.e2=pmet;
  PD.f1=alca;
  PD.g=Lsov;
%  PD.r=Perdu;
%  PD.s=Perdb;
  PD.t=Tefv;
 elseif idyn==2
  PD.g=Lsov;
 end
global n_pa n_ve  
  ngra=[n_pa n_ve];
  PD.nefgra=ngra;
%' PD', keyboard

ierr=0;
if iriga==1
 if npfr~=0
  modaut1
 else
  fris=fr;
 end
 if exist('gg0')==0
  ierr=1;
 end
% disp(' dopo modaut')
% pausak
 if ierr==0
  if ifnm==0
   riga
  else
   save mop
   mop1
   disp(' verifica qui ')
   keyboard
%   autdiff1
  end
   if ifeed==1
    subrife
   else
    Pfi=0;
    Pfi1=0;
    pi0=0;
   end
 else
  dnup=0;
  ae=0;
  dnupv=0;
  aev=0;
  gsoav=0;
  fsoav=0;
  asoav=0;
  Nth=0;
  fTras=0;
  Pfi=0;
  Pfi1=0;
  pi0=0;
  Ww=0;
  Is=0;
  Nno=0;
  Psi=0;
  Alfp=0;
  fTras=0;
  nsp=0;
 end
else
 dnup=0;
 ae=0;
 dnupv=0;
 aev=0;
end
if ifp>=-2
 pausak
end
if ierr==1
 disp(' ierr=1 ')
 pausak
end

if icampi>0
 if exist('Gsov')
 if length(Gsov)>0
  gsov=Gsov;
  fsov=Fsov;
  ADo=ADout;
  nazim=Nazim;
  nrad=Nrad;
  M_EDc=EDc;
  M_EDo=EDout;
 end
 end
end

if istop_fie==1
 ifp=ifpsalvo;
end 
 if isfield(Ps,'ifpstop')==1
  if Ps.ifpstop==1
   ifp=ifpsalvo;   
  end
 end  

%if ispeed==1
%  eval(['!rd ',Dsav]);
%end
%' fine subolt: controllo ', keyboard
