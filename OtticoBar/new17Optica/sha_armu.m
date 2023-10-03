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
Delta=-.3;
%Delta=0;
cce=[0 4 j*2 3+j*2];

npfia=101;
fi1=linspace(1e-7,pi-1e-7,npfia)';
fi=[fi1; fi1+pi];

fid=linspace(0,2*pi,npfia)';
fi=fid(1:end-1);


cces=cce;
ccesav=cce;

fice=find(real(cce)>0 & imag(cce)>0);
fice0=find(abs(cce)~=0);

cce0=ccesav(fice).';
cce=ccesav(fice).';
cce=[cce; -cce];
cce=[cce; conj(cce)];

if length(fice0)>0
 ficed=find((real(cce)==0) & (abs(real(cce))+abs(imag(cce))~=0));
 if length(ficed)>0
  cce=[cce; ccesav(ficed) ];
 end
 ficed=find((imag(cce)==0) & (abs(real(cce))+abs(imag(cce))~=0));
 if length(ficed)>0
  cce=[cce; ccesav(ficed) ];
 end
end

fice00=find((abs(cce)==0));
 if length(fice00)>0
  cce=[cce; ccesav(fice00) ];
 end


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
 for nh=1:length(cce)
  ro=cce(nh)+R*(1+Delta*cos(4*fi)).*exp(j*fi);
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
for nh=1:length(cce0)
 ro=cce0(nh)+R*(1+Delta*cos(4*fi)).*exp(j*fi);
 M=abs(ro);
 F=angle(ro);
 fiF=find(F<0);
 F(fiF)=F(fiF)+2*pi;
 Rv=[Rv ro];
 Mv=[Mv M];
 Fv=[Fv F];
end
Fv=unwrap(Fv);

for nh=1:length(cce0)
 Fd=Fv(:,nh);
 if length(find(Fd>0))<length(Fd)/2
  Fv(:,nh)=Fd+2*pi;
 end
end

dro=1e-4;
sR=size(Mv);
rdu=reshape(Mv,1,prod(sR));
srdu=sort(rdu);
fisrdu=find(diff([0 srdu])>=dro);
rod=srdu(fisrdu);
ro=(rod(1:end-1)+rod(2:end))/2;
ro=[rod(1)+1e-6 ro rod(end)-1e-6];
%figure, plot(ro,'.')


rax=max(ro)*1.05;
ng=length(ro);
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

%  ro=rv'/asav*avero;
  rv=ro'/avero*asav;


Rlo=ones(ng,6*length(cce0))*NaN;

%if ifp>=0
% fiN=figure;
%end
fiq=linspace(0,pi/2,51);
cir=exp(j*fiq);
nsol=zeros(size(rv));

icdisp=100;
coic=0;
for ic=1:length(cce0)
%for ic=3
 mafi=1;
 fiinv0=ones(ng,4)*1.5;
 rop=cce0(ic)+R*(1+Delta*cos(4*fi)).*exp(j*fi);

 roo=abs(rop);
 [maro,ima]=max(roo);
 [miro,imi]=min(roo);
 punacc=find(ro+dro<maro & ro-dro>miro);

  if ic==icdisp
   hj=figure;
   hf=figure;
  end

% puni=find(abs(ro-miro)<=dro);
% pup=nsol(puni);
% nsol(puni)=pup+2;
% Rlo(puni,(1:2)+pup)=rop(imi);
% puni=find(abs(ro-maro)<=dro);
% pup=nsol(puni);
% nsol(puni)=pup+2;
% Rlo(puni,(1:2)+pup)=rop(ima);

 for ipcel=punacc
% for ipcel=25:30
   rvi=ro(ipcel);
   fii=find(rvi==roo);

   y=roo-rvi;
   y1=y([end 1:end-1]);
   y2=y(1:end);
   yp=y1.*y2;
   fiia1=find(yp<0);
   fiiad=[fii' fiia1'];
   lfi=length(fiiad);
%   if lfi>4
%    lfi
%    keyboard
%   end
%   if lfi==3
   if fix(lfi/2)-lfi/2~=0
    lfi
    keyboard
    [fa,ifa]=sort(angle(rop(fiiat)));
    dufi=ifa([1 3]);
    fiiat=fiiad(dufi);
   else
    fiiat=fiiad;
   end
   if ic==icdisp
    fiiat
    pausak
   end

   roz=[];
   fiia=[];
   for kf=1:length(fiiat)
    fii=fiiat(kf);
    if rvi~=roo(fii)
     fiia=[fiia fii];
     if fii~=1
      pu=fii+[-1 0];
     else
      pu=[length(rop) 1];
     end
     an=angle(rop(pu));
     ri=roo(pu);
     cfip=polyfit(ri,an,length(ri)-1);
     ri0=roo(pu);
     dx=max(abs(diff(ri0)));
     xfip=linspace(min(ri0)-dx/2,max(ri0)+dx/2,20);
     yfip=polyval(cfip,xfip);
     fip=polyval(cfip,rvi);
     if ic==icdisp
      figure(hf); plot(ri,an,'*',xfip,yfip,rvi,fip,'ro'), pausak
     end
     roz=[roz rvi*exp(j*fip)];
    else
     roz=[roz rop(fii)];
    end
   end
   lz=length(roz);
   pup=nsol(ipcel);
%   if pup>0
%    pup
%    keyboard
%   end
   Rlo(ipcel,(1:lz)+pup)=roz;
   nsol(ipcel)=pup+lz;
   if ic==icdisp
    figure(hj); plot(rop,'.'), axis equal, hold on
    plot(cir*rvi,'c'),  plot(rop([fiiat]),'ro'),
    plot(roz,'m.'),
    pausak, clf
   end

 end
% figure, plot(Rlo,'.'), axis equal,  pausak

end

Rsa=Rlo;
%Rlo=[Rlo -Rlo];
%Rlo=[Rlo conj(Rlo)];

 [Flo,ipu]=sort(angle(Rlo),2);
 for k=1:length(ipu)
  Mlo(k,:)=abs(Rlo(k,ipu(k,:)));
 end
 figure, polar(Flo,Mlo,'.'), pausak
 mans=max(nsol);
 pui=[1:2:mans];
 puf=[2:2:mans];
 Floi=Flo(:,pui);
 Flof=Flo(:,puf);
 Mloi=Mlo(:,pui);
 Mlof=Mlo(:,puf);
 figure, plot(Mloi,Floi,'g.',Mlof,Flof,'r.'),
 pausak


fasid=Floi;
fasud=Flof;
Mid=Mloi;
Mud=Mlof;

fasid=[fasid pi-Flof];
fasud=[fasud pi-Floi];

fasid=[fasid pi+Floi];
fasud=[fasud pi+Flof];

fasid=[fasid 2*pi-Flof];
fasud=[fasud 2*pi-Floi];

for kf=1:3
 Mid=[Mid Mloi];
 Mud=[Mud Mlof];
end


 figure, plot(Mid,fasid,'g.',Mud,fasud,'r.'),
 pausak


r=ro;
ru=Mv;
xt=Fv;


if ifp>1
 figure, plot(r,fasid,'g.',r,fasud,'r.',ru,xt),
 pausak
end

muv=[0:2:2*(nubesu)];
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
figure, plot(ro,AB)

cce=cces;
igint=igintsa;
