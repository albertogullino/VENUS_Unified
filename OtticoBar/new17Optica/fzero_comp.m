%'fzero_comp', keyboard

ifpsal=ifp;

ifp=ifceck_zero;
if ifceck_zero==1
 ifp=-10;
end


igrad=1;

if igrad==1
 nr=2;
% nr=5;
else
 nr=3;
% nr=5;
end 


%[fpe,fpm]=f_mulut(0);
% tolto provvisorio
%[fpe,fpm]=f_mulut(0);

icomp_grat=0;
icomp_grat=1;
if  iscan==1
ng=11;
nl=13;
z0=j*nime*100;
str=.005;
%str=str0;
nl=13;
nl=21;
sti=nime*100;
save sal
lav=real(z0)+linspace(-str,str,nl);
Lam0=lambda_cen+lime;
%lav=linspace(-str,str,nl);
niv=imag(z0)+linspace(-sti,sti,ng);
niv=imag(z0)+linspace(0,sti,ng);
%niv=0;

lavm=repmat(lav',1,ng);
nivm=repmat(niv,nl,1);
z=lavm+j*nivm;

z=lav+j*nime;

ifun=0;




'dopo promap'
keyboard
%keyboard
%keyboard
icomp_grat=1;
 [fpe,fpm]=f_mulut(z);
%  figure, plot(z,abs(fpe),z,abs(fpm)) 
 if itetmt==1
  fp=fpe;
 else
  fp=fpm;
end
'dopo mulut', keyboard
% fp=f_mul(z);

%[prev,idu1]=min(abs(fp));
%[pre,idu2]=min(prev);
%z0v=z0;
%z0=z(idu1(idu2),idu2);
%fp0=fp(idu1(idu2),idu2);

[pre,idu1]=min(abs(fp));
z0v=z0;
z0=z(idu1);
fp0=pre;

str=diff(lav(1:2))*5;
str=diff(lav(1:2));
%str=str0;

%ipfi1=idu1+[-1 0 1];
%fpf=fp(ipfi1);
%zf=z(ipfi1);
%cof=polyfit(zf,fpf,1);
%z00=roots(cof);

%sti=diff(niv(1:2));

if ifp==-10

figure, semilogy(lav,abs(fp),real(z0),pre,'wo'), pausak
if iscan==1
 iok=input(' OK ? 0,no ');
 if length(iok)==0
  iok=1;
 end 
 if iok==0
  ag=ginput(1);
  z0=ag(1);
  str=str0/20;
 end
end
end

precf=1e-4;
precx=1e-7;
precg=1e-3;
% ifp=-10
imold=1;
if imold==0
 nl=3;
 ng=3;
else
 nl=5;
 ng=5;
end
zerev=[];
pre=100;
icit=0;
tic
if ifp==-10
 hh=figure;
end
% icomp_grat=1;
'qui', keyboard
while pre>precf
icit=icit+1;
if str<1e-3
% icomp_grat=0;
end
if icit>2
% icomp_grat=0;
% icomp_grat=0;
end
if imold==1
lav=real(z0)+linspace(-str,str,nl);
niv=imag(z0)+linspace(-sti,sti,ng);

lavm=repmat(lav',1,ng);
nivm=repmat(niv,nl,1);
z=lavm+j*nivm;

else
  %'qui', keyboard
   if min(size(fp))>1
    refp=mean(real(fp)');
    isiz=2;
   else
    refp=real(fp);
    isiz=1;
   end
   [du,imi]=min(abs(refp-lav));
   co=polyfit(lav,refp,1);
   zere=roots(co);
   [du,imi]=sort(abs(refp-lav));
  if isiz==2
   imfp=mean(imag(fp(imi(1:2),:)));
   co=polyfit(niv,imfp,1);
   zeim=roots(co);
  else
   imfp=imag(fp);
   zeim=nime;
  end


  if length(zerev)==1
   str=abs(zerev-zere);
   sti=abs(zeimv-zeim);
   zerev=zere;
   zeimv=zeim;
  else
   str=str/2;
   sti=sti/2;
   zerev=zere;
   zeimv=zeim;  
  end
  lav=zere+linspace(-str,str,nl);
  niv=zeim+linspace(-sti,sti,ng);
% 'zere', zere, zeim

  lavm=repmat(lav',1,ng);
  nivm=repmat(niv,nl,1);
  z=lavm+j*nivm;
end

ifun=0;
 [fpe,fpm]=f_mulut(z);
 if itetmt==1
  fp=fpe;
 else
  fp=fpm;
end
'dopo mulut', keyboard
[prev,idu1]=min(abs(fp));
[pre,idu2]=min(prev);
z0v=z0;

if imold==1
z0=z(idu1(idu2),idu2);

%pausak
fp0=abs(fp(idu1(idu2),idu2));
if ifp==-10
z0
fp0
end
str=diff(lav(1:2));
sti=diff(niv(1:2));
else
z0=zere+j*zeim;
pre=min(min(abs(fp)));
end
zeim=imag(z0);
if ifp==-10
 figure(hh), semilogy(lav,abs(fp),real(z0),pre,'wo'), 
 sti/zeim, 
 str,
 pre
pausak
end

 if sti/zeim<precg | str<precx
  icit
  break
 end
end
toc
if ifp==-10
figure, semilogy(lav,abs(fp),real(z0),pre,'wo'), 
'residuo', abs(fp0)
pausak
end


else  %scan

load contzer
iczer=iczer+1;
save contzer icscan iczer
ifun=itetmt;
%' nime ', keyboard
zr0=nime;
zv_0=5e-5*linspace(-1,1,nr);
%zv_0=1e-4*linspace(-1,1,nr);
%zv_0=1e-3*linspace(-1,1,nr);
if igrad==2
 zvr=zr0*[.5 1 1.5];
else
 zvr=zr0*[.2  1.];
end
clear zi zr
for ked=1:length(zvr)
zv=zv_0+j*zvr(ked);
 [fpe,fpm]=f_mulut(zv);
 fp_conf=mean(abs(fpe));
 if itetmt==1
  fpro=fpe;
 else
  fpro=fpm;
 end

fr=real(fpro);
fi=imag(fpro);
cor=polyfit(zv_0,fr,igrad);
coi=polyfit(zv_0,fi,igrad);
if igrad==2
 zrdu=roots(cor);
 [du,is]=min(abs(zrdu));
 zr(ked)=zrdu(is);
 zidu=roots(coi);
 [du,is]=min(abs(zidu));
 zi(ked)=zidu(is);
else
 zrdu=roots(cor);
 zr(ked)=zrdu;
 zidu=roots(coi);
 zi(ked)=zidu;
end
 if ifp==-10
  figure, plot(zv,real(fpro),zr(ked),0,'wo',zv,imag(fpro),zi(ked),0,'wo'),
  title(' dlv_0 ')
  pausak
 end 
end  %ked

cor=polyfit(zvr,(zr-zi),igrad);
corf=polyfit(zvr,zr,igrad);

if igrad==2
 zugd=roots(cor);
 [du,is]=min(abs(zugd-mean(zvr)));
 zug=zugd(is);
else
 zugd=roots(cor);
 zug=zugd;
end

fug=polyval(corf,zug);
 if ifp==-10
  figure, plot(zvr,zr,'.-',zvr,zi,'.-',zug,fug,'wo'),
  title(' dlv_0 (nI) ')
  pausak
end

zver=fug+j*zug;
[fpe,fpm]=f_mulut(zver);
 if itetmt==1
  fpro0=abs(fpe);
 else
  fpro0=abs(fpm);
 end
iwcont=0;
%' prec', keyboard
precf=.01;
if lsce==1
precf=.001;

end
while fpro0>precf & iwcont<2
zv_0=fug+1e-5*linspace(-1,1,nr);
if igrad==2
 zvr=zug*[.9 1 1.1];
else
 zvr=zug*[.8+iwcont*.08  1.2-iwcont*.08];
end
iwcont=iwcont+1;
icomp_grat=1;
clear zi zr
for ked=1:length(zvr)
zv=zv_0+j*zvr(ked);
 [fpe,fpm]=f_mulut(zv);
 if itetmt==1
  fpro=fpe;
 else
  fpro=fpm;
 end
%fpro=f_mulast(zv);
fr=real(fpro);
cor=polyfit(zv_0,fr,igrad);
zrdu=roots(cor);
if igrad==2
 [du,is]=min(abs(zrdu-fug));
 zr(ked)=zrdu(is);
else 
 zr(ked)=zrdu;
end
fi=imag(fpro);
coi=polyfit(zv_0,fi,igrad);
zidu=roots(coi);
if igrad==2
 [du,is]=min(abs(zidu-fug));
 zi(ked)=zidu(is);
else
 zi(ked)=zidu;
end
 if ifp==-10
  figure, plot(zv,real(fpro),zr(ked),0,'wo',zv,imag(fpro),zi(ked),0,'wo'),
  title(' seconda iterazione ')
  pausak
 end 

 
end  %ked
cor=polyfit(zvr,(zr-zi),igrad);
corf=polyfit(zvr,zr,igrad);
%zug=roots(cor);
zugd=roots(cor);
if igrad==2
 [du,is]=min(abs(zugd-mean(zvr)));
 zug=zugd(is);
else
 zug=zugd;
end


fug=polyval(corf,zug);

zver=fug+j*zug;
%fver=f_mulast(zver)




%' veri refi', keyboard
 if ifp==-10
  figure, plot(zvr,zr,'.-',zvr,zi,'.-',zug,fug,'wo'),
  title(' seconda iterazione, Guadagno')

  pausak
 end 
zver=fug+j*zug;
[fpe,fpm]=f_mulut(zver);
 if itetmt==1
  fpro0=abs(fpe);
 else
  fpro0=abs(fpm);
 end
 iwcont 
 fpro0
% keyboard
 if fpro0<precf
  break
 end 
end % if fpro
z0=zver;
end %iscan

ifp=ifpsal;