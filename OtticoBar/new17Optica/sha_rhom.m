%ifp=10
%if istrumix==1
% aS=a;
% bS=b;
%
% a=aloc;
% b=bloc;
%end
%igint=0;
%'entro sha_rhom', keyboard
if ~exist('iBWnew')
 iBWnew=0;
end

if igraef_new>=0 & iBWnew>=1
' strato con variazione trasversale e reticolo'
  ainf=0;
  asup=aloc;

   ab=asup;
   aa=ainf;
   ng=fix(ab*max(KK)*200);
  if igint==1
   [rv,wi]=gauleg(aa,ab,ng);
  else
   rv=linspace(aa,ab,ng);
  end
  sgimp=1;
  muv=[0:pasnu:2*nubesu];
  fafa=pi/2;
  AB=zeros(size(muv));
  AB(1)=fafa;
  funzAB
  
  r=rv/(kcav);
  ifalso=-1;   % non devo fare parte circolare analitica!
  % gia fatto
  igf=0
  if igf==0
  fireBW=find(shavet==-8);
  par_BW=radii.array{fireBW(1)}{13};
  DC=par_BW.DC;
  if ifp==-10
  ' qui anisot', keyboard
  end
  
  if DC>0
   g_n1=par_BW.n1;
   g_n2=par_BW.n2;
   er1=g_n1.^2;
   er2=g_n2.^2;
   e_ve=1./(DC/er1+(1-DC)/er2);
   e_pa=DC*er1+(1-DC)*er2;
   Period=par_BW.period*1000;
   lambda1000=lambda*1000;
   thickness=Litot(fireBW-1)*1000;
   if fireBW<fiQW
    nin=nitot(fireBW,1);
    nout=nitot(fireBW-2,1);
   else
    nout=nitot(fireBW,1);
    nin=nitot(fireBW-2,1);  
   end
   NModes=11;
  
   [neff,nBW,flagt] = f_EffectiveIndex(Period,lambda1000,DC,nin,nout,g_n1,g_n2,thickness,NModes);
   if iBWnew==1
    e_pa=nBW(1)^2;
    e_ve=nBW(2)^2;
   else
    e_pa=neff(1)^2;
    e_ve=neff(2)^2;  
   end
%  n_ve=sqrt(e_ve);
%  n_pa=sqrt(e_pa);  
%  if ifp==-10
%   'qui ALB', keyboard
  else  % anisotropia Liquid Crystal
%  ' qui anisot LC', keyboard  
    nany=nitot(fireBW-1,:);
    neff(1)=nany(1);
    neff(2)=par_BW.period;
%    neff(2)=neff(1);
%    'NO ANY CONFINATA!!!!!!!',    pausak
   
    e_pa=neff(1)^2;
    e_ve=neff(2)^2;

  end
%  end 
end %igf
%  neff(2)=neff(1);
  enne_med=sqrt((e_ve+e_pa)/2);
  Sig_eps=(e_ve+e_pa)/(rr^2)-2;
  Sig_eps=0;
  Del_eps=(e_ve-e_pa)/(rr^2);
  Del_eps=1;
  %Sig_eps=(e_ve+e_pa)/(2*rr^2)-1;
  %Del_eps=(e_ve-e_pa)/(2*rr^2);
  
%             depepfz=-rr^2/riv(ianu)^2+rr^2/rext^2;  
  Delz_eps=-rr^2/e_pa; 
    if ifp==-10
     ' qui anisot LC', keyboard 
    end
 ioold=0;
 if ioold==1
  if ipolar==-1
  Sig_eps=2*e_pa/rr^2-2;
  Delz_eps=-rr^2/e_pa;
  else 
  Sig_eps=2*e_ve/rr^2-2;
  Delz_eps=-rr^2/e_ve;  
  end
  Del_eps=0;
 end


  alons=alon;
  
%  r=ru*R00;
%  'qui RETICOLO', keyboard

 return
end

alimm=max(alim);

imeshae=0; %=1 e` simile a ellisse (come area considerata), =0 mantiene area costante
if istrumix==0
 if abs(b-a)>10
  ng=fix(abs(b-a)*alimm*20);
 else
  ng=fix(b*alimm*3*3);
 end
else
 if abs(bloc-aloc)>10
  ng=fix(abs(bloc-aloc)*alimm*20);
 else
  ng=fix(bloc*alimm*20);
 end
end
if ng<50
 ng=50;
end
disp(' shaoc')
%keyboard

nv_sh=4;
fia=2*pi/nv_sh;
nfia=501;
xt=linspace(0,2*pi,nfia);

if istrumix==0
 rapax=bvero/avero;
 averodi=avero;
elseif istrumix==1
 rapax=bloc/aloc;
 averodi=aloc/kcav;
 avero=aloc/kcav;
 bvero=bloc/kcav;
end
%dae=(rapax-1)/2;
%'sha_rom ', keyboard

rapax=bvero/avero;
dae=-(1-rapax)/(1+rapax)*(1+dap);
R00=bvero/(1+dae+dap);
%ru=R00*(1+dap*cos(nv_sh*xt)+dae*cos(2*xt));

if imeshae==0
 ru=1+dap*cos(nv_sh*xt)+dae*cos(2*xt);
else
 ru=1+dap*cos(nv_sh*xt)+dae*cos(2*xt)+dae;
 ru=1+dap*cos(nv_sh*xt)+dae*cos(2*xt)+dae;
end
 ruc=1+0*cos(nv_sh*xt);

if ifp>-3
%if ifp==-10
 figure, polar(xt,ru*R00,'r'),
 hold on, polar(xt,ruc*R00,'g'), title(' sha-oxi '), pausak
 figure, plot(xt,ru-1,'r'), pausak
end

% Rvd=R00*ru.*exp(j*xt);
% figure, plot(Rvd), grid, axis equal
%if ifp==-10
% figure, polar(xt,ru*R00,'r'),
% hold on, polar(xt,ruc*R00,'g'), title(' sha-rhom '),
% drawnow
%end

%a=min(R00*ru)*kcav;

%asav=R00*kcav;
asav=min(R00*ru)*kcav;
alon=asav;

rmi=min(ru)*R00*(1+1e-5);
%rmi=min(ru);
[du,fiad]=min(abs(adis-a));
%if abs(adis-a*rmi)>1e-6
% adis(fiad)=a*rmi;
% a=a*rmi;
%end

rma=max(ru)*R00;
%rma=max(ru);
%rmima=[rmi rma];
%sgsha=1;
%[du,idu]=min(abs(rmima-ru(1)));
%if idu==1
% sgsha=-1;
%end

  ainf=rmi*asav;
  asup=rma*asav;

  ainf=rmi*kcav;
  asup=rma*kcav;

%rmi=min(ru*R00);
%rma=max(ru*R00);

%ainf=rmi*kcav0;
%asup=rma*kcav0;
%asav=ainf;

%  if sgsha==1
   ab=asup;
   aa=ainf;
%  else
%   aa=asup;
%   ab=ainf;
%  end
  if igint==1
   [rv,wi]=gauleg(aa,ab,ng);
  else
   rv=linspace(aa,ab,ng);
  end
  r=rv/(kcav);
%  r=ru*R00;
%'cont rv', keyboard
  sgimp=ones(size(r));
%  fis=find(r<1);
%  sgimp(fis)=-1;

cx(1)=8*dap;
cx(2)=2*dae-8*dap;
cx(1:2)=cx(1:2)*R00;
fas=ones(length(r),6)*NaN;
for k=1:length(r)

%  cx(3)=(1+dap-dae)-r(k);
  cx(3)=(1+dap-dae)*R00-r(k);
%' harom', keyboard

 rd=roots(cx);
 fi=find(imag(rd)==0 & real(rd)>=0 & real(rd)<=1);
 if length(fi)~=0
  cf=sqrt(rd(fi));
  fas(k,1:length(fi))=acos(cf)';
 end
end

fast=[fas pi-fas fas+pi 2*pi-fas];
s=size(fast);
fasi=ones(length(r),16)*NaN;
fasu=ones(length(r),16)*NaN;
for k=1:length(r)
  fd=fast(k,:);
  fin=find(isnan(fd)==0);
  lf=length(fin);
  if lf~=0
   fdv=sort(fd(fin));
   fii=1:2:lf;
   fiu=2:2:lf;
   fasi(k,1:length(fii))=fdv(fii);
   fasu(k,1:length(fiu))=fdv(fiu);
  end
end
ruR00=ru*R00;
if ruR00(1)>rmi
%if ru(1)>rmi
 [du,fim]=min(fasi(:,1));
 faa=0*fasu(:,1);
 ps=fim+1:length(faa);
% if du>1
  faa(ps)=NaN*faa(ps);
% end
 fasid=[faa fasu];
 faa=2*pi*fasi(:,1)./fasi(:,1);
% if du>1
  faa(ps)=NaN*faa(ps);
% end
 fasud=[fasi faa];
 fasdu=fasid;
 fasid(ps,:)=fasud(ps,:);
 fasud(ps,:)=fasdu(ps,:);
else
 fasid=fasi;
 fasud=fasu;
end
if ifp>1
 figure, plot(r,fasid,'g.',r,fasud,'r.',ruR00,xt),
 pausak
end
muv=[0:pasnu:2*nubesu];
%muv=[0:2:18];
AB=ones(length(r),length(muv))*NaN;
CD=AB;
im=0;
for mu=muv
im=im+1;
 for k=1:length(r)
  fid=fasid(k,:);
  fud=fasud(k,:);
  fini=find(isnan(fid)==0);
  finu=find(isnan(fud)==0);
  if length(fini)~=0 & length(fini)==length(finu)
   if mu==0
    Ad=sum(fud(finu)-fid(fini));
    Cd=0;
   else
    Ad=sum(sin(mu*fud(finu))-sin(mu*fid(fini)))/mu;
    Cd=-sum(cos(mu*fud(finu))-cos(mu*fid(fini)))/mu;
   end
%   pausak
   AB(k,im)=Ad/2;
   CD(k,im)=Cd/2;
  end
 end
 fiA=find(abs(AB(:,im))<1e-7);
 if length(fiA)>ng/2
  AB(:,im)=AB(:,im)*0;
 end
end
if ifp>1
% figure, plot(r,A(:,1:6)), hold on, plot(r,A(:,7:length(muv)),'.-'),
% title(' y0, m2, c4, r6, g8, b10, y12, m14, c16, r18, g20, b22')
 if length(muv)>6
  figure, plot(r,AB(:,1:6)), hold on, plot(r,AB(:,7:length(muv)),'.-'),
 else
  figure, plot(r,AB),
 end
 title(' y0, m2, c4, r6, g8, b10, y12, m14, c16, r18, g20, b22')
 pausak
end

  if ifp>10
   figure
  end

 for imu=pimu
  jmu=imu-meun;
  mu=mbv(imu);
  for inu=pimu
   jnu=inu-meun;
   nu=mbv(inu);
%   if (nu+mu)/2-fix((nu+mu)/2)==0
    dmn=abs(mu-nu);
    sdmn=sign(mu-nu);
%    if sdmn==0
%     keyboard
%    end
    fim=find(dmn==muv);
    if length(fim)==1
     mfatd=AB(:,fim);
     A(:,jmu,jnu)=AB(:,fim);
     if ipolar==0
      C(:,jmu,jnu)=sdmn*CD(:,fim);
     end
    else
     disp('errore A mu in sha_oxi ')
     pausak
    end
    dmn=abs(mu+nu);
    fim=find(dmn==muv);
    if length(fim)==1
     mfats=AB(:,fim);
     B(:,jmu,jnu)=AB(:,fim);
     if ipolar==0
      D(:,jmu,jnu)=CD(:,fim);
     end
    else
     disp('errore B mu in sha_oxi ')
     pausak
    end
   end  %if
%  end
 end
%' e questa ', keyboard  
  
  
return

'vecchio caso simmetrico'
clear A B
 for imu=pimu
  jmu=imu-meun;
  mu=mbv(imu);
  for inu=pimu
   jnu=inu-meun;
   nu=mbv(inu);
   if (nu+mu)/2-fix((nu+mu)/2)==0
    dmn=abs(mu-nu);
    fim=find(dmn==muv);
    if length(fim)==1
     mfatd=AB(:,fim);
     A(:,jmu,jnu)=AB(:,fim);
    else
     disp('errore A mu in sha_oxi ')
     pausak
    end
    dmn=abs(mu+nu);
    fim=find(dmn==muv);
    if length(fim)==1
     mfats=AB(:,fim);
     B(:,jmu,jnu)=AB(:,fim);
    else
     disp('errore B mu in sha_oxi ')
     pausak
    end
   end  %if
  end
 end

%a=asav;

%if istrumix==1
% a=aS;
% b=bS;
%end

if istrumix==1
 aloc=ainf;
 bloc=asup;
end

%rv=r*kcav0;
%' e questa ', keyboard
