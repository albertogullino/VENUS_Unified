ifp_shav=ifp;
%ifp=4
fice=find(real(cce)>0 | imag(cce)>0);
fice0=find(abs(cce)==0);
ficet=find(real(cce)>=0 & imag(cce)>=0);
ficet=1;
fi=linspace(0,2*pi,4*npfia+1)';

 Rvd=[];
 for nh=1:length(cce)
  ro=sha_gen(Pshi,fi,nh);
  Rvd=[Rvd; ro];
 end
if ifp>=1
%if ifp>=1 | ifp==-10

% fip=find(real(Rvd)>=0 & imag(Rvd)>=0);
% Rac=Rvd(fip);
% Rvt=[Rac; -Rac];
% Rvt=[Rvt; conj(Rvt)];
 Rvt=Rvd;
 figure, plot(Rvt,'.'), hold on,
 plot(Rvd,'r.'), axis equal,
 grid
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
 ro=sort([ro; Ma+1e-6; Mi-1e-6 ]);
end

if ifp>2
figure, plot(ro,'.'), pausak
%keyboard
end

rax=max(ro)*1.05;
ng=length(ro);
axe=[-rax rax -rax rax];

if ng<20
 ng=20;
end

asav=a;

%  ro=rv'/asav*avero;
%  rv=ro'/avero*asav;
  rv=ro';
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
fiq=linspace(0,pi*2,51);
cir=exp(j*fiq);
nsol=zeros(size(rv));

icdisp=500;
%icdisp=1;
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

% fiok=find(real(ropd)>=-1e-12 & imag(ropd)>=-1e-12);
% rop=ropd(fiok);
 rop=ropd;

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
% for ipcel=10
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
%     ' UNO '
    elseif real(rop(fii))==0
     roz=[roz j*rvi];
%     ' DUE '
    else
%     ' TRE '
     if rvi~=roo(fii)
      fiia=[fiia fii];
      if fii~=1
       pu=fii+[-1 0];
      else
       pu=[length(rop) 1];
      end
      an=unwrap(angle(rop(pu)));
%      an=(angle(rop(pu)));
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

%   if length(find(real(roz)<-1e-10))>0
%    '<0'
%    keyboard
%   end

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
if ifp>=0
 figure, polar(Flo,Mlo/kcav0,'.'), pausak
end
 mans=max(nsol);
 pui=[1:2:mans];
 puf=[2:2:mans];
 Floi=Flo(:,pui);
 Flof=Flo(:,puf);
 Mloi=Mlo(:,pui);
 Mlof=Mlo(:,puf);
% figure, plot(Mloi,Floi,'g.',Mlof,Flof,'r.'),
% pausak


r=ro;

fad=Flo;
fi=find(Flo<0);
fad(fi)=fad(fi)+2*pi;
fad=sort(fad')';
fi=find(fad==0);
fad(fi)=NaN;
%figure, plot(r,fad,'.')
 Floi=fad(:,puf);
 Flof=fad(:,pui);
%figure, plot(r,Floi,'g.',r,Flof,'ro')

rvet=repmat(r,1,min(size(fad)));

fadv=reshape(fad,prod(size(fad)),1);
fadr=reshape(rvet,prod(size(fad)),1);
[Fadv,is]=sort(fadv);
Fadr=fadr(is);
iva=find(isnan(Fadv)==0);
Fadr=Fadr(iva);
Fadv=Fadv(iva);

fi=find(diff(Fadr)==0);
[du,ima]=max(Fadv);
pii=[fi(end)+1:ima];
Fadv(pii)=Fadv(pii)-2*pi;
[Fadv,is]=sort(Fadv);
Fadr=Fadr(is);


fi=find(diff(Fadr)==0);
%fi=[1; fi; ima];
fi=[fi; ima];
% figure, plot(Fadr,Fadv,'-',Fadr(fi),Fadv(fi),'r.')
% pausak

[mi,im]=min(Fadr(fi));
fisa=fi;
sufa=ones(size(fisa));
%if im>1
 fi=[fisa(im:end); fisa(1:im-1)];
 sufa=[sufa(im+1:end); sufa(1:im)*0];
%end

fipro=[fi+1 fi];
[du,imfi]=max(fipro);
fipro(imfi)=1;
fii=fipro(:,1);
fiu=[fipro(2:end,2); fipro(1,2)];

% figure,
% for ipl=1:8
%  plot(Fadr,Fadv,'-',Fadr(fii(ipl)),Fadv(fii(ipl)),'g.',...
%                            Fadr(fiu(ipl)),Fadv(fiu(ipl)),'r.')
%  pausak
% end

%fii=[1; sort([fi(1:2:end-1); fi(2:2:end-1)+1]) ];
%fiu=[fi(end-1)+1; fi(end); sort([fi(1:2:end-2)+1; fi(2:2:end-2)]) ];


% figure, plot(Fadr,Fadv,'-',...
% Fadr(fii),Fadv(fii),'g.',Fadr(fiu),Fadv(fiu),'ro')
% pausak


%pii=[];
%piu=[];
%for ii=1:2:length(fiu)-1
% pii=[pii fii(ii):fii(ii+1)];
% piu=[piu fiu(ii):fiu(ii+1)];
%end

  pii=[];
  piu=[];
  fashii=[];
  fashiu=[];
  for ii=1:2:length(fiu)-1
   pii=[pii fii(ii):fiu(ii)];
   fashii=[fashii ones(size(fii(ii):fiu(ii)))*sufa(ii)];
   piu=[piu fii(ii+1):fiu(ii+1)];
   fashiu=[fashiu ones(size(fii(ii+1):fiu(ii+1)))*sufa(ii+1)];
  end
% figure, plot(Fadr,Fadv,'-',Fadr(pii),Fadv(pii)-fashii'*2*pi,...
%                'g.',Fadr(piu),Fadv(piu)-fashiu'*2*pi,'r.')
% pausak

Fadv(pii)=Fadv(pii)-fashii'*2*pi;
Fadv(piu)=Fadv(piu)-fashiu'*2*pi;

Fri=Fadr(pii);
Ffi=Fadv(pii);

Fru=Fadr(piu);
Ffu=Fadv(piu);

clear fasid fasud
fasid=ones(length(r),mans)*NaN;
fasud=fasid;


for kfa=1:length(r)
 fiid=find(Fri==r(kfa));
 fiud=find(Fru==r(kfa));
 fasid(kfa,1:length(fiid))=Ffi(fiid)';
 fasud(kfa,1:length(fiud))=Ffu(fiud)';
end


%if diff(fasid(1:2,1))<0
% fad=fasid;
% fasid=fasud;
% fasud=fad+2*pi;
%end
if mean(fasud(:,1))<0
 fasid=fasid+2*pi;
 fasud=fasud+2*pi;
end

%if ifp<=-10 | ifp>0
if ifp>=0
 figure, plot(r,(fasid),'g.',r,(fasud),'r.'),
 pausak
end
ifp=ifp_shav;