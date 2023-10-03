%if istrumix==1
% aS=a;
% bS=b;
%
% a=aloc;
% b=bloc;
%end

imeshae=0; %=1 e` simile a ellisse (come area considerata), =0 mantiene area costante
if istrumix==0
 if abs(b-a)>10
  ng=fix(abs(b-a)*alim*20);
 else
  ng=fix(b*alim*3*3);
 end
else
 if abs(bloc-aloc)>10
  ng=fix(abs(bloc-aloc)*alim*20);
 else
  ng=fix(bloc*alim*3*3);
 end
end
if ng<20
 ng=20;
end
%disp(' shaoc')
%keyboard

nv=4;
fia=2*pi/nv;
nfia=501;
xt=linspace(0,2*pi,nfia);

if istrumix==0
 rapax=bvero/avero;
 averodi=avero;
elseif istrumix==1
 rapax=bloc/aloc;
 averodi=aloc/kcav;
end
%dae=-0.03;
dae=(rapax-1)/2;
%dap=0.05;
%dap=0.0;
if imeshae==0
 ru=1+dap*cos(nv*xt)+dae*cos(2*xt);
else
 ru=1+dap*cos(nv*xt)+dae*cos(2*xt)+dae;
 ru=1+dap*cos(nv*xt)+dae*cos(2*xt)+dae;
end
 ruc=1+0*cos(nv*xt);
if ifp>-3
 figure, polar(xt,ru*averodi,'r'),
 hold on, polar(xt,ruc*averodi,'g'), title(' sha-oxi '), pausak
 figure, plot(xt,ru-1,'r'), pausak
end
if ifp==-10
 figure, polar(xt,ru*averodi,'r'),
 hold on, polar(xt,ruc*averodi,'g'), title(' sha-oc '),
 drawnow
end

asav=a;
rmi=min(ru);
[du,fiad]=min(abs(adis-a));
if abs(adis-a*rmi)>1e-6
 adis(fiad)=a*rmi;
 a=a*rmi;
end

rma=max(ru);
%rmima=[rmi rma];
%sgsha=1;
%[du,idu]=min(abs(rmima-ru(1)));
%if idu==1
% sgsha=-1;
%end


  ainf=rmi*asav;
  asup=rma*asav;
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
  r=rv/asav;
  sgimp=ones(size(r));
%  fis=find(r<1);
%  sgimp(fis)=-1;

cx(1)=8*dap;
cx(2)=2*dae-8*dap;
fas=ones(length(r),6)*NaN;
for k=1:length(r)
 if imeshae==0
  cx(3)=1-r(k)+dap-dae;
 else
  cx(3)=1-r(k)+dap;
 end
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
if ru(1)>rmi
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
 figure, plot(r,fasid,'g.',r,fasud,'r.',ru,xt),
 pausak
end
muv=[0:2:2*nubesu];
%muv=[0:2:18];
AB=ones(length(r),length(muv))*NaN;
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
   else
    Ad=sum(sin(mu*fud(finu))-sin(mu*fid(fini)))/mu;
   end
%   pausak
   AB(k,im)=Ad/2;
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
 aloc=a;
 bloc=bloc*rmi;
end
