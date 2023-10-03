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
Delta=-.5;
Delta=0;
cce=[j*2+2  4+j*7 7];


npfia=61;
np00=31;

fi=linspace(0,2*pi,npfia)';
fi00=linspace(0,2*pi,2*npfia)';

fi1q=find(fi<pi/2);
fi1q0=find(fi00<pi/2);
fi4q=find(fi>3*pi/2);

cces=cce;
ccesav=cce;

fice=find(real(cce)>0 | imag(cce)>0);

cce0=ccesav(fice).';

fice0=find(abs(cce)==0);


pice=find(real(cce)>0 & imag(cce)>0);
pice0=find(abs(cce)~=0);

ccep=ccesav(pice).';
ccep=[ccep; -ccep];
ccep=[ccep; conj(ccep)];

if length(pice0)>0
 ficed=find((real(cce)==0) & (abs(cce)~=0));
 if length(ficed)>0
  ccep=[ccep; cce(ficed); -cce(ficed)];
 end
 ficed=find((imag(cce)==0) & (abs(cce)~=0));
 if length(ficed)>0
  ccep=[ccep; cce(ficed); -cce(ficed)];
 end
end


%if ifp>=1
%if ifp>=1 | ifp==-10
 Rvd=[];
 for nh=1:length(ccep)
  ro=ccep(nh)+R*(1+Delta*cos(4*fi)).*exp(j*fi);
  Rvd=[Rvd ro];
 end
 figure, polar(angle(Rvd),abs(Rvd))
% figure, plot(real(Rvd),imag(Rvd)), axis equal,
 drawnow
 pausak
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

if length(fice0)>0
 rod0=abs(cce(fice0)+R*(1+Delta*cos(4*fi00)).*exp(j*fi00));
else
 rod0=[];
end

dro=1e-4;
sR=size(Mv);
rdu=reshape(Mv,1,prod(sR));
srdu=sort([rdu rod0']);
fisrdu=find(diff([0 srdu])>=dro);
rod=srdu(fisrdu);
ro=(rod(1:end-1)+rod(2:end))/2;

if length(fice0)>0
 ro=sort([ro max(Mv)+1e-6 min(Mv)-1e-6 max(rod0)+1e-6]);
 ro00=linspace(0,min(rod0),np00);
 ro=[ro00(1:end-1) ro(1:end)];
else
 ro=sort([ro max(Mv)+1e-6 min(Mv)-1e-6 ]);
end

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


Rlo=ones(ng,6*length(cce))*NaN;

%if ifp>=0
% fiN=figure;
%end
fiq=linspace(0,pi/2,51);
cir=exp(j*fiq);
nsol=zeros(size(rv));

icdisp=100;
coic=0;
for ic=1:length(cce)
%for ic=3

 if abs(cce(ic))~=0
  ropd=cce(ic)+R*(1+Delta*cos(4*fi)).*exp(j*fi);
 else
  ropd=cce(ic)+R*(1+Delta*cos(4*fi00)).*exp(j*fi00);
 end
 fiok=find(real(ropd)>=0 & imag(ropd)>=0);
 if length(fiok)<length(ropd)
  if imag(cce(ic))==0 & real(cce(ic))~=0
   fi0=[pi 0];
   ropdui=cce(ic)+R*(1+Delta*cos(4*fi0)).*exp(j*fi0);
   radd=linspace(ropdui(1),ropdui(2),10).';
   rop=[ropd(fiok); radd];
  elseif imag(cce(ic))~=0  & real(cce(ic))==0
   fi0=[pi/2 3*pi/2];
   ropdui=cce(ic)+R*(1+Delta*cos(4*fi0)).*exp(j*fi0);
   radd=linspace(ropdui(1),ropdui(2),10).';
   rop=[ropd(fi1q); radd; ropd(fi4q)];
  else  %centro in 0
%   fi0=[pi/2];
%   ropdui=cce(ic)+R*(1+Delta*cos(4*fi0)).*exp(j*fi0);
%   radd=linspace(ropdui(1),0,10).';
%   rop=[ropd(fi1q0); radd];
%
%   fi0=[0];
%   ropdui=cce(ic)+R*(1+Delta*cos(4*fi0)).*exp(j*fi0);
%   radd=linspace(0,ropdui(1),10).';
%   rop=[rop; radd];

   rop=ropd(fi1q0);
  end
 else
  rop=ropd(fiok);
 end


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
% for ipcel=25:35
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
%   if ic==icdisp
%    fiiat
%    pausak
%   end

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
   if length(find(real(roz)<-1e-10))>0
    '<0'
    keyboard
   end
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

if length(fice0)>0
 fir=find(ro<=min(rod0));
 fasid(fir,1)=0;
 fasud(fir,1)=2*pi;
 Mid(fir,1)=ro(fir)';
 Mud(fir,1)=ro(fir)';
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