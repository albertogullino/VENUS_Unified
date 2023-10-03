igint=0;
R=Rx;


%fi0=[];
%r0a=[];
%fi0add=0;
%fi1=pi/2*[0  2];
%r1a=2.5*ones(size(fi1));
%r1a=3.5*ones(size(fi1));
%fi0=[fi0 fi1];
%r0a=[r0a r1a];

cces=cce;
icirc0=0;
fic=find(abs(cce)==0);
if length(fic)>0
 fic=find(abs(cce)~=0);
 cce=cces(fic);
 icirc0=1;
end

%if ifp>=1
if ifp>=1 | ifp==-10
 Rvd=[];
 for nh=1:length(cces)
  ro=cces(nh)+R*exp(j*fi);
  Rvd=[Rvd ro];
 end
% figure, plot(real(Rvd),imag(Rvd)), axis equal,
% drawnow
% if ifp>=1, pausak, end
end


r0a=abs(cce);
fi0=angle(cce);

Mv=[];
Fv=[];
Rv=[];
for nh=1:length(cce)
 ro=cce(nh)+R*exp(j*fi);
 M=abs(ro);
 F=angle(ro);
 fiF=find(F<0);
 F(fiF)=F(fiF)+2*pi;
 Rv=[Rv ro];
 Mv=[Mv M];
 Fv=[Fv F];
end
Fv=unwrap(Fv);

for nh=1:length(cce)
 Fd=Fv(:,nh);
 if length(find(Fd>0))<length(Fd)/2
  Fv(:,nh)=Fd+2*pi;
 end
end

%Fvi=[];
%Fvu=[];
%Mvi=[];
%Mvu=[];
%ifm=find(fi0>pi);
%fi00=fi0;
%fi00(ifm)=fi00(ifm)-2*pi;
%for nh=1:length(fi0)
% Fd=Fv(:,nh);
% Md=Mv(:,nh);
% Fdm=fi00(nh);
% ifi=find(Fd<Fdm);
% Fvi=[Fvi Fd(ifi)];
% Mvi=[Mvi Md(ifi)];
% ifu=find(Fd>Fdm);
% Fvu=[Fvu Fd(ifu)];
% Mvu=[Mvu Md(ifu)];
%end


rax=max(r0a)+max(R);
ng=fix(a*alim*8*rax/avero);
if ng<20
 ng=20;
end
if ifp>=1
 figure, plot(real(Rv),imag(Rv)), axis square, axis([-rax rax -rax rax]), pausak
 figure, plot(Fv,Mv),
 axis([min(min(Fv)) max(max(Fv)) 0 rax]), pausak
 figure, plot(Mv,Fv), pausak
end

asav=a;
rmi=(min(r0a)-R)/avero;
rma=(max(r0a)+R)/avero;

ainf=rmi*asav;
asup=rma*asav;
ab=asup;
aa=ainf;

%  if igint==1
%   [rv,wi]=gauleg(aa,ab,ng);
%  else
   rv=linspace(aa,ab,ng);
%  end

  ro=rv'/asav*avero;
  wi=[0 diff(rv)];


fivf=ones(ng,length(cce))*NaN;
fivi=fivf;
for ic=1:length(cce)
 r0ai=r0a(ic);
 fiinv0=(ro.^2-(R^2+r0ai^2))/(2*r0ai*R);
 fiv=find(fiinv0<=1 & fiinv0>=-1);
 fiinv=fiinv0(fiv);
 lefi=[1:length(fiv)];
 lefu=[1:length(fiv)]+length(fiv);
 fiinv=[fiinv; flipud(fiinv)];
 ofiinv=[-ones(size(fiv)); ones(size(fiv))];
 ro0=ro(fiv);
% ic
% keyboard

 fi0i=0;
 andu=fi0i+acos(fiinv).*ofiinv;
 arg=(r0ai*sin(fi0i)+R*sin(andu))./(r0ai*cos(fi0i)+R*cos(andu));
 ata=unwrap(2*atan(arg))/2;

% figure, subplot(311), plot(andu),
% subplot(312), plot(arg),
% subplot(313), plot(ata), pausak
% rov=[rov ro0];
 fi0i=fi0(ic);
 fivf(fiv,ic)=flipud(ata(lefu))+fi0i;
 fivi(fiv,ic)=ata(lefi)+fi0i;
end

if ifp>=1
figure, plot(ro,fivi), hold on, plot(ro,fivf,'--')
pausak
end

fasid=fivi;
fasud=fivf;

r=ro;
ru=Mv;
xt=Fv;


if ifp>1
 figure, plot(r,fasid,'g.',r,fasud,'r.',ru,xt),
 pausak
end

muv=[0:2:2*nubesu];
%muv=[0:2:18];
AB=zeros(length(r),length(muv));
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
%   if abs(Cd)>1e-6
%    Cd
%    pausak
%   end
   AB(k,im)=Ad/2;
   CD(k,im)=Cd/2;
  end
 end
 fiA=find(abs(AB(:,im))<1e-14);
 AB(fiA,im)=AB(fiA,im)*0;
 fiA=find(abs(CD(:,im))<1e-14);
 CD(fiA,im)=CD(fiA,im)*0;

% fiA=find(abs(AB(:,im))<1e-7);
% if length(fiA)>ng/2
%  AB(:,im)=AB(:,im)*0;
% end
% fiA=find(abs(CD(:,im))<1e-7);
% if length(fiA)>ng/2
%  CD(:,im)=CD(:,im)*0;
% end
end
if ifp>1
 if length(muv)>6
  figure, plot(r,AB(:,1:6)), hold on, plot(r,AB(:,7:length(muv)),'.-'),
 else
  figure, plot(r,AB),
 end
 title(' AB: y0, m2, c4, r6, g8, b10, y12, m14, c16, r18, g20, b22')
 pausak
 if length(muv)>6
  figure, plot(r,CD(:,1:6)), hold on, plot(r,CD(:,7:length(muv)),'.-'),
 else
  figure, plot(r,CD),
 end
 title(' CD: y0, m2, c4, r6, g8, b10, y12, m14, c16, r18, g20, b22')
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

cce=cces;
