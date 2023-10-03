disp(' sha_arj')
igintsa=igint;
igint=0;
%keyboard
ifalso=-1

sgimp=1;
if istrumix==1
 avero=aloc/kcav0;
end

R=avero;
R=1.6
%disp(' R'), keyboard
Delta=dap;
Delta=-.2;
%Delta=0;

npfia=301;
fi1=linspace(1e-7,pi-1e-7,npfia)';
fi=[fi1; fi1+pi];

fid=linspace(0,2*pi,npfia)';
fi=fid(1:end-1);


cces=cce;
icirc0=0;
fic=find(abs(cce)==0);
if length(fic)>0
 fic=find(abs(cce)~=0);
 cce=cces(fic);
 icirc0=1;
end

%if ifp>=1
%if ifp>=1 | ifp==-10
 Rvd=[];
 for nh=1:length(cces)
  ro=cces(nh)+R*(1+Delta*cos(4*fi)).*exp(j*fi);
  Rvd=[Rvd ro];
 end
 figure, polar(angle(Rvd),abs(Rvd))
% figure, plot(real(Rvd),imag(Rvd)), axis equal,
 drawnow
 if ifp>=1, pausak, end
%end


r0a=abs(cce);
fi0=angle(cce);

Mv=[];
Fv=[];
Rv=[];
for nh=1:length(cce)
 ro=cce(nh)+R*(1+Delta*cos(4*fi)).*exp(j*fi);
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


rax=max(r0a)+max(R)+2*abs(Delta);
ng=fix(a*alim*8*rax/avero);
ng=fix(101*rax/avero);
axe=[-rax rax -rax rax];
if ng<20
 ng=20;
end
if ifp>=1
 figure, plot(real(Rv),imag(Rv)), axis square, axis(axe), pausak
 figure, plot(Fv,Mv),
 axis([min(min(Fv)) max(max(Fv)) 0 rax]), pausak
 figure, plot(Mv,Fv), pausak
end

asav=a;
rmi=(min(r0a)-R-2*abs(Delta))/avero;
rma=(max(r0a)+R+2*abs(Delta))/avero;

ainf=rmi*asav;
asup=rma*asav;
ab=asup;
aa=ainf;
  if igint==1
   [rv,wi]=gauleg(aa,ab,ng);
  else
   rv=linspace(aa,ab,ng);
  end

  ro=rv'/asav*avero;



fivf=ones(ng,4*length(cce))*NaN;
fivi=fivf;

%if ifp>=0
% fiN=figure;
%end

coic=0;
for ic=1:length(cce)
%for ic=1
 mafi=1;
 fiinv0=ones(ng,4)*1.5;
 rop=cce(ic)+R*(1+Delta*cos(4*fi)).*exp(j*fi);

    rop1=abs(cce(ic))+R*(1+Delta*cos(4*fi)).*exp(j*fi);
 [maro,ima]=max(abs(rop1));
 [miro,imi]=min(abs(rop1));
 vmima(ic,1)=miro;
 vmima(ic,2)=maro;
 ropmi=rop(imi);
 ropma=rop(ima);
% if ifp>=0
%  figure(fiN), plot(rop,'y'), hold on, plot(ropmi,'w*'), plot(ropma,'r*'),
%  axis(axe), axis equal,
% end
 r0ai=r0a(ic);
 if Delta==0
  fiinv0=(ro.^2-(R^2+r0ai^2))/(2*r0ai*R);
 else
  cte(1)=(8*R*Delta)^2;
  cte(2)=0;
  cte(3)=-128*(R*Delta)^2;
  cte(4)=16*R*r0ai*Delta;
  cte(5)=64*(Delta*R)^2+16*(1+Delta)*Delta*R^2;
  cte(6)=-16*R*r0ai*Delta;
  cte(7)=-16*R^2*Delta*(1+Delta);
  cte(8)=2*r0ai*R*(1+Delta);

  firo=find(ro>=miro & ro<=maro);
  for kr=firo'
   cte(9)=-(ro(kr).^2-(R^2*(1+Delta)^2+r0ai^2));
   rote=roots(cte);
%   pausak
   fii=find(imag(rote)==0);
   fiac=find(abs(rote(fii))<=1);
   if length(fiac)>4
    ro(kr)
    an=acos(rote(fii(fiac)))
    rop1=abs(cce(ic))+R*(1+Delta*cos(4*fi)).*exp(j*fi);
    cir=ro(kr)*exp(j*fi);
    figure, plot(rop1,'y'), hold on, plot(cir,'r'),
    axis equal, grid,
    pausak
   end
   lfia=length(fiac);
   if lfia>0
    fiinv0(kr,1:length(fiac))=rote(fii(fiac))';
    if lfia>mafi
     mafi=lfia;
    end
   end
  end
 end

% verifica poli
% cte(9)=(R^2*(1+Delta)^2+r0ai^2);
% rop=r0ai+R*(1+Delta*cos(4*fi)).*exp(j*fi);
% fx=fi;
% x=cos(fx);
% y=polyval(cte,x);
% figure, plot(fx,sqrt(y),fi,abs(rop)), pausak

 for ks=1:mafi
  coic=coic+1;
  fiv=find(abs(fiinv0(:,ks))<=1);
  fiinv=fiinv0(fiv,ks);

  lefi=[1:length(fiv)];
  lefu=[1:length(fiv)]+length(fiv);
  fiinv=[fiinv; flipud(fiinv)];
  ofiinv=[-ones(size(fiv)); ones(size(fiv))];

%  fi0i=0;
  fi0i=fi0(ic);
  andu=fi0i+acos(fiinv).*ofiinv;
  Randu=R*(1+Delta*cos(4*andu));
  arg=(r0ai*sin(fi0i)+Randu.*sin(andu))./(r0ai*cos(fi0i)+Randu.*cos(andu));
  ata=unwrap(2*atan(arg))/2;

%  fi0i=fi0(ic);
  fi0i=0;
  fivf(fiv,coic)=flipud(ata(lefu))+fi0i;
  fivi(fiv,coic)=ata(lefi)+fi0i;
  figure, plot(flipud(ata(lefu)),'.'), hold on, plot(ata(lefi),'r.')
  pausak
 end  %ks

end

% if ifp>=0
%  pausak
% end

%fiN=find(isnan(fivi(:,mafi))==0);
%pui=[1:2:2*mafi];
%puf=[2:2:2*mafi];
%for kr=fiN';
%  dui=fivi(kr,:);
%  duf=fivf(kr,:);
%  dut=[dui duf];
%  fiNt=find(isnan(dut)==0);
%  duta=sort(dut(fiNt));
%  fivi(kr,1:mafi)=duta(pui);
%  fivf(kr,1:mafi)=duta(puf);
%end

for kr=1:length(ro);
  dui=fivi(kr,:);
  duf=fivf(kr,:);
  dut=[dui duf];
  fiNt=find(isnan(dut)==0);
  lfiNt=length(fiNt);
  if lfiNt>0
   pui=[1:2:lfiNt];
   puf=[2:2:lfiNt];
   duta=sort(dut(fiNt));
   fivi(kr,1:lfiNt/2)=duta(pui);
   fivf(kr,1:lfiNt/2)=duta(puf);
  end
end

if ifp>=0
figure, plot(ro,fivi,'.g'), hold on, plot(ro,fivf,'.r')
pausak

% rp=ro;
% Fi=fivi;
% Ff=fivf;
% save sa rp Fi Ff

% load sa
% figure, plot(ro,fivi), hold on, plot(ro,fivf,'--')
% hold on, plot(rp,Fi,'.'), hold on, plot(rp,Ff,'.')
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
igint=igintsa;
