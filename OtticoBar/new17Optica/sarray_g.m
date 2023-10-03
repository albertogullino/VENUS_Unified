
fice=find(real(cce)>0 | imag(cce)>0);
fice0=find(abs(cce)==0);
ficet=find(real(cce)>=0 & imag(cce)>=0);


%if ifp>=1
if ifp>=1 | ifp==-10
 Rvd=[];
 for nh=1:length(cce)
  ro=sha_gen(Pshi,fi,nh);
%  keyboard
  Rvd=[Rvd; ro];
 end

 fip=find(real(Rvd)>=0 & imag(Rvd)>=0);
 Rac=Rvd(fip);
 Rvt=[Rac; -Rac];
 Rvt=[Rvt; conj(Rvt)];
 figure, plot(Rvt,'.'), hold on,
 plot(Rvd,'r.'), axis equal,
 if ifp>=0, pausak, end
 hold off
end

Mv=[];
Ma=[];
Mi=[];
sfic=size(fice);
if sfic(1)>sfic(2)
 fice=fice';
end
%for nh=fice
for nh=ficet
 ro=sha_gen(Pshi,fi,nh);
 M=abs(ro);
 Ma=[Ma; max(M)];
 Mi=[Mi; min(M)];
 Mv=[Mv; M];
end

if length(fice0)>0
 rod0=abs(sha_gen(Pshi,fi00,fice0));
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

asav=a;

%  ro=rv'/asav*avero;
  rv=ro'/avero*asav;
  wi=[0 diff(rv)];
%  if igint==1
%   [rv,wi]=gauleg(aa,ab,ng);
%  else
%   rv=linspace(aa,ab,ng);
%  end


Rlo=ones(ng,6*length(cce))*NaN;

%if ifp>=0
% fiN=figure;
%end
fiq=linspace(0,pi/2,51);
cir=exp(j*fiq);
nsol=zeros(size(rv));

icdisp=500;
sfic=size(ficet);
if sfic(1)>sfic(2)
 ficet=ficet';
end

for ic=ficet

 if abs(cce(ic))~=0
  ropd=sha_gen(Pshi,fi,ic);
 else
  ropd=sha_gen(Pshi,fi00,ic);
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
% figure, plot(Mloi,Floi,'g.',Mlof,Flof,'r.'),
% pausak


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

% figure, plot(Mid,fasid,'g.',Mud,fasud,'r.'),
% pausak


r=ro;

if ifp>1
 figure, plot(r,fasid,'g.',r,fasud,'r.'),
 pausak
end
