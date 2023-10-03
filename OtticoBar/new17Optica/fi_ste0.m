colordef('black')
clear
clear global
close all
 co(1:2,1)='w-';
 co(1:2,2)='m-';
 co(1:2,3)='c-';
 co(1:2,4)='g-';
 co(1:2,5)='r-';
 co(1:2,6)='b-';
 co(1:3,7)='y--';
 co(1:3,8)='m--';
 co(1:3,9)='c--';
 co(1:3,10)='g--';
 co(1:3,11)='r--';
 co(1:3,12)='b--';

 co1(1:2,1)='w.';
 co1(1:2,2)='m.';
 co1(1:2,3)='c.';
 co1(1:2,4)='g.';
 co1(1:2,5)='r.';
 co1(1:2,6)='b.';
%
r_pil=2.5;

N=0;
T=0;

if T~=0 | N~=0
load tempdistrib
disp('Temp')
keyboard
load tempdistrib
%load carrier2mA
delr=18/100;
delz=7.656276829/100;
[R,Z] = meshgrid(delr:delr:18,7.6563:-delz:delz);
%[C,h] = contour(R,Z,TempDistrib,10);
%clabel(C,h)
%colormap(cool)
dndt=2.3e-4;
nmed=3.3;
Tp=2*nmed*dndt*(TempDistrib-297);
ro=R(1,:);
zeta=Z(:,1);
sro=12;
Np=exp(-(ro/sro).^2);
%Np=spline(ra,N,ro*1e-6);
else
 ro=0;
 zeta=0;
end

if N~=0
 N=Np;
end

if T~=0
 T=Tp;
end


%ird=[ 1:16:100];
%figure, plot(ro,T(ird,:)'), xlabel(' ro ')
%figure, plot(zeta,T(:,ird)), xlabel(' zeta ')
%figure, plot(ro,dndt*T(ird,:)'), xlabel(' ro ')
%figure, plot(ro,N), xlabel(' ro ')
%figure, plot(zeta,dndt*T(:,ird)), xlabel(' zeta ')
%pausak


shape0='sha_oc';
i2D=3;
sha=1;
mmvet=[0 1];
sha=4;
iLP=0;
mmvet=[3 2]+1;
numodiacc=2;     % numero di modi da accoppiare (sopra, sotto)

mmvet=[3 2];
mmvet=[3];
numodiacc=1;     % numero di modi da accoppiare (sopra, sotto)
mmvet=[3];
numodiacc=1;     % numero di modi da accoppiare (sopra, sotto)

%mmvet=[8 9];
%numodiacc=5;     % numero di modi da accoppiare (sopra, sotto)
nomes='dp';

%iLP=1;
%sha=1;
%mmvet=[0:14];
%numodiacc=0;     % numero di modi da accoppiare (sopra, sotto)

iany=0;  % 1: anisotropia planare, 2: confinata (effetto elettro-ottico )
%iany=0;  % 1: anisotropia planare, 2: confinata (effetto elettro-ottico )
ianys=0;   % strain

if iLP==1
 ianys=0;
 iany=0;
end

ifp=0;
ifp=-5;
%ifp=-10;
%ifp=-3;

Frisi0=0.5e-3;
Frisu0=0.9e-3;
Ndisp0=5;
nk1max=20;
%nk1max=30;
alim=0.20;
idyn=0;

Frisi0=0.4e-3;
Frisu0=0.8e-3;

Frisi0=0.4e-3;
Frisu0=0.8e-3;
Ndisp0=3;
nk1max=20;
alim=0.2;
numodiacc=2;

%nk1max=10;
%alim=0.15;
%numodiacc=1;
ifp=-4;

nvar='s_stein';

%% ifp:  -10 alcuni stop campi
%% ifp:  -5  no display campi e no stop
%% ifp:  -4  display campi e no stop

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% calcolo e plot modi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iraff=1;
iraff=-2;

if iraff==0
 disp(' per ora iraff = 0 ')
end

t0=clock;
%xroI1=linspace(0,3,30);
%xroI2=linspace(3,6,10);
%xroI=[xroI1 xroI2(2:length(xroI2))];

ploma=2*r_pil;
Np=60;
npfi0=15;

%%%%%%%%
 ipurstep=0;  % =0 step costante

% r=linspace(0,ploma,Np);

global Dapv ras igamveb igamveu GGbext GGuext
igamveb=1;   %=1, gamma vero,  =0, riflessione GGext
igamveu=1;
GGbext=1;
%GGbext=1;
GGuext=1;

nv=4;
fia=2*pi/nv;
nfia=4*(npfi0)+1;
xt=linspace(0,2*pi,nfia);
global Dapv ras
rapax=0.85;
ras=rapax;
if sha==1
 rapax=1;
end
dae=(rapax-1)/2;
dap=0.04;
Dapv=dap;
ru=(1+dap*cos(4*xt)+dae*cos(2*xt));
 r=linspace(0,ploma,Np-1);

% ihol=.5;
% pj_ulm
% close all

 xroI=r;
 xro=r;

%percm=r_pil/ploma;
%
%% grigliato non omogeneo
%
%x=linspace(0,5,Np-1);
%y=-exp(-.5*(x-7*percm*.9).^2)+.2*exp(atan(x));
%ym=min(y);
%S=y-ym+.05;
%Sm=sum(S);
%r=[0 cumsum(S)/Sm*ploma];
%figure, plot(diff(r),'.')
%figure, plot(r,'.'),
%pausak

%xro=r;







lfi_inp=2*npfi0;
fimaxi=pi;
iplan=1;

fileeps='Gres_jo2';
Dla=-.005;


disp(' call caloptg '),
nell=3.559124817121504e+000 -6.740211839941768e-004i;

par_inv=[4.5/2 2.25/2 nell;
         4.5/2 2.25/2 1;
         1.8/2 7.2/2 1];

%par_inv=[];
sp=size(par_inv);
nloo=max(sp(1),1);

for loop=1:nloo
 if sp(1)>0
  par_in=par_inv(loop,:);
 else
  par_in=[];
 end

 [Eqw,Eout,xro,fian,lambda,delf,gain,ord,nrAzim,Cu,PaDy]=...
 calopt_m(N,T,ro,zeta,xroI,fimaxi,lfi_inp,par_in,r_pil,nvar,...
 shape0,mmvet,numodiacc,Frisi0,Frisu0,alim,Ndisp0,nk1max,...
 iplan,iraff,idyn,i2D,iLP,sha,iany,ianys,ifp,Dla,fileeps,nomes);



 tim=etime(clock,t0);
 ore=fix(tim/3600);
 minuti=fix((tim-ore*3600)/60);
 secondi=(tim-ore*3600-minuti*60);
 format short
 disp(' elapsed time (hours, minutes, seconds) = ');
 disp([ ore minuti secondi]);

 lamm(:,loop)=lambda./(1+delf);
 demm(:,loop)=delf;
 pamm(:,loop)=ones(size(delf))*loop;
 gamm(:,loop)=gain';
 tymm(:,loop)=PaDy.l';
end
fiX=find(tymm>0);
fiY=find(tymm<0);
coGH=3e5/lambda;
demg=demm-min(min(demm));
figure, plot(demg(fiX)*coGH,gamm(fiX),'rx'), hold on,
plot(demg(fiY)*coGH,gamm(fiY),'go'), grid,
title('red x: xpol       green o: ypol')
Bir=(demm(fiY)-demm(fiX))*coGH;
Dic=(gamm(fiY)-gamm(fiX));
tab(:,1)=pamm(fiX);
tab(:,2)=Bir;
tab(:,3)=Dic;
tab
bir=tab(:,2);
dic=tab(:,3);
figure, plot(bir,dic,'ro')
xlabel(' fy-fx')
ylabel(' Gy-Gx')
famo=(bir/min(bir)-1)/10;
famol=1+max(famo);
for im=1:length(bir);
   gmod=dic(im);
   lmod=bir(im);
   sym=[fchar(im)];
   ht=text(lmod*famol,gmod,sym); set(ht,'fontsize',8);
end %im


%stat
keyboard
if ifp>-4
 close all
 keyboard
end


figure, plot(lamm,gain,'r.')
hold on
famo=(lamm/min(lamm)-1)/10;
famol=1+max(famo);
for im=1:length(lamm);
   gmod=gain(im);
   lmod=lamm(im);
   inua=ord(1,im);
   inur=ord(2,im);
   sym=[fchar(inua) '-' fchar(inur)];
   ht=text(lmod*famol,gmod,sym); set(ht,'fontsize',8);
end %im

% pcam

if i2D==3
 xcu=Cu.x;
 ycu=Cu.y;
 zcu=Cu.z;
 lg=length(gain);
 nsi=floor(sqrt(lg));
 nsu=ceil(sqrt(lg));
 if nsi*nsu<lg
  nsi=nsu;
 end
 XP=xro'*cos(fian);
 YP=xro'*sin(fian);
 aax=1.2*max([max(xcu) max(ycu)]);
 tyPmod=PaDy.l;
 figure
 pograp=[ 10 50 1250 900];  set(gcf,'Position',pograp);
 for k=1:lg
    if tyPmod(k)>0
     Emoto=reshape(Eout.x(:,:,k),length(xro),length(fian));
    else
     Emoto=reshape(Eout.y(:,:,k),length(xro),length(fian));
    end
    subplot(nsi,nsu,k),
    surf(XP,YP,Emoto),
    shading('interp'),
    view(0,90),
    axis square, axis equal, axis off, axis([-1 1 -1 1]*aax),
     stri=[' gain = ',fcharc(gain(k),7)];
     text(-aax,1.4*aax,stri);
     stri=[' freq = ',fcharc(delf(k),7)];
     text(-aax,1.2*aax,stri);
     if tyPmod(k)>0
      stri='X_{pol}';
     else
      stri='Y_{pol}';
     end
     text(1.2*aax,0,stri);
    hold on, fill3(xcu,ycu,zcu,'w')
 end

end

idiffu=1;

if iLP==1
 nomefi='deb_ulmLP';
else
 nomefi='deb_ulm';
end

eval(['save ' nomefi]);

disp(' sto entrando in dynvet ')
if ifp>-4
 keyboard
end


