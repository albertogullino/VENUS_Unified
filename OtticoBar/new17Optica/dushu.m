igintsa=igint;
igint=0;
%keyboard
ifalso=-1

sgimp=1;
if istrumix==1
 avero=aloc/kcav0;
end

isamecell=1;

if isamecell==1
 ncex=3;
 ncey=3;
 dxce=1;
 dyce=1;
 ncell=ncex*ncey;
 Rx=2;
 Ry=2;
 dcex=2*Rx+dxce;
 dcey=2*Ry+dyce;
 Del=.3;
 sh_ty=1;

 Rvetx=repmat(Rx,1,ncell);
 Rvety=repmat(Ry,1,ncell);
 Delta=repmat(Del,1,ncell);
 sh_type=repmat(sh_ty,1,ncell);

 fcex=([1:ncex]-ncex/2-.5)*dcex;
 fcey=([1:ncey]-ncey/2-.5)*dcey;
 Fcx=ones(size(fcey'))*fcex;
 Fcy=fcey'*ones(size(fcex));
 cced=Fcx+j*Fcy;
 cce=reshape(cced,1,ncell);
 figure, plot(cce,'o'), pausak
% if is_even(ncex)==1 & is_even(ncey)==1
% elseif is_even(ncex)==1 & is_even(ncey)==0
% elseif is_even(ncex)==0 & is_even(ncey)==1
% elseif is_even(ncex)==0 & is_even(ncey)==0
% end

else

 cce=[0 j*3+2  4+j*7 7 7*j];
 Rvetx=[.4  1      .5   2  .8];
 Rvety=[.2  2      .5   3  1];
 Delta=[ 0 -.2     .5  .8 -.05];
 sh_type=[1 2       1   1   2];

% Delta=[ 0 .2     .5  .8  .05];
% sh_type=[1 1       1   1   1];

end

Psh{1}=sh_type;
Psh{2}=Rvetx;
Psh{3}=Rvety;
Psh{4}=Delta;
Psh{5}=cce;



npfia=11*numodiacc;
np00=31;

fi=linspace(0,2*pi,4*npfia+1)';
fi00=linspace(0,2*pi,8*npfia+1)';


fice=find(real(cce)>0 | imag(cce)>0);
fice0=find(abs(cce)==0);
ficet=find(real(cce)>=0 & imag(cce)>=0);


%if ifp>=1
%if ifp>=1 | ifp==-10
 Rvd=[];
 for nh=1:length(cce)
  ro=sha_gen(Psh,fi,nh);
  Rvd=[Rvd; ro];
 end
 figure, plot(Rvd,'.'), axis equal
 drawnow
 pausak
 if ifp>=1, pausak, end
%end

 fip=find(real(Rvd)>=0 & imag(Rvd)>=0);
 Rac=Rvd(fip);
 Rvt=[Rac; -Rac];
 Rvt=[Rvt; conj(Rvt)];
 figure, plot(Rvt,'.'), axis equal
 pausak


Mv=[];
Ma=[];
Mi=[];
for nh=fice
 ro=sha_gen(Psh,fi,nh);
 M=abs(ro);
 Ma=[Ma; max(M)];
 Mi=[Mi; min(M)];
 Mv=[Mv; M];
end

if length(fice0)>0
 rod0=abs(sha_gen(Psh,fi00,fice0));
else
 rod0=[];
end

dro=1e-12;
srdu=sort([Mv; rod0]);
fisrdu=find(diff([0; srdu])>=dro);
rod=srdu(fisrdu);
ro=(rod(1:end-1)+rod(2:end))/2;

if length(fice0)>0
 ro=sort([ro; Ma+1e-6; Mi-1e-6; max(rod0)+1e-6]);
 ro00=linspace(0,min(rod0),np00)';
 ro=[ro00(1:end-1); ro(1:end)];
else
 ro=sort([ro; Ma+1e-6; Mi-1e-6 ]);
end

%figure, plot(ro,'.')
%keyboard


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

icdisp=500;

for ic=ficet

 if abs(cce(ic))~=0
  ropd=sha_gen(Psh,fi,ic);
 else
  ropd=sha_gen(Psh,fi00,ic);
 end

 fiok=find(real(ropd)>=-1e-12 & imag(ropd)>=-1e-12);
 rop=ropd(fiok);

  if imag(cce(ic))==0 & real(cce(ic))~=0
   setmet=1;
  elseif imag(cce(ic))~=0  & real(cce(ic))==0
   setmet=2;
  else
   setmet=0;
  end


 if setmet==2
%   rop=rop([end/2+2:end 1:end/2+1]);
   rop=rop([end/2+1:end 1:end/2]);
 end

 roo=abs(rop);

 [maro,ima]=max(roo);
 [miro,imi]=min(roo);
 punacc=find(ro+dro<maro & ro-dro>miro)';

  if ic==icdisp
   hj=figure;
   hf=figure;
  end


 for ipcel=punacc
   rvi=ro(ipcel);
   fii=find(rvi==roo);

   y=roo-rvi;
   if setmet==0
    y1=y([end 1:end-1]);
    y2=y(1:end);
    sup=0;
   elseif setmet>0
    y1=y([1:end-1]);
    y2=y(2:end);
    sup=1;
   end

   yp=y1.*y2;
   fiia1=find(yp<0)+sup;
   fiiad=[fii' fiia1'];
   lfi=length(fiiad);
   if lfi==1
    if setmet==1
     fiiat=[fiiad length(rop)];
    elseif setmet==2
     fiiat=[fiiad 1];
    end
   else
    fiiat=fiiad;
   end

   roz=[];
   fiia=[];
   for kf=1:length(fiiat)
    fii=fiiat(kf);
    if imag(rop(fii))==0
     roz=[roz rvi];
    elseif real(rop(fii))==0
     roz=[roz j*rvi];
    else
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
   end
   lz=length(roz);
   pup=nsol(ipcel);

   if length(find(real(roz)<-1e-10))>0
    '<0'
    keyboard
   end

   Rlo(ipcel,(1:lz)+pup)=roz;
   nsol(ipcel)=pup+lz;
   if ic==icdisp
    figure(hj); plot(rop,'.'), axis equal, hold on
    plot(cir*rvi,'c'),
    roplr=real(rop(fiiat));
    ropli=imag(rop(fiiat));
    rozplr=real(roz);
    rozpli=imag(roz);
    plot(roplr,ropli,'ro'),
    plot(rozplr,rozpli,'m.'),
    pausak,
    clf
   end

 end
% figure, plot(Rlo,'.'), axis equal,  pausak

end

Rsa=Rlo;

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
 Mid(fir,1)=ro(fir);
 Mud(fir,1)=ro(fir);
end

for k=1:3
 fiN=find(abs(fasid-pi/2*k)<=1e-12);
 fasid(fiN)=NaN;
% fiN=find(fasud==pi/2*k);
 fiN=find(abs(fasud-pi/2*k)<=1e-12);
 fasud(fiN)=NaN;
end

 figure, plot(Mid,fasid,'g.',Mud,fasud,'r.'),
 pausak


r=ro;

if ifp>1
 figure, plot(r,fasid,'g.',r,fasud,'r.'),
 pausak
end

muv=[0:2:2*(nubesu)];
muv=[0:2:30];
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
 mfiA=mean(abs(AB(:,im)))
 if mfiA<2e-3
  AB(:,im)=0;
  ABac(im)=0;
 else
  ABac(im)=1;
 end
% fiA=find(abs(CD(:,im))<1e-14);
% CD(fiA,im)=CD(fiA,im)*0;
 mfiA=mean(abs(CD(:,im)));
 if mfiA<1e-4
  CD(:,im)=0;
 end

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
figure, plot(ro,AB,ro,CD)

cce=cces;
igint=igintsa;
