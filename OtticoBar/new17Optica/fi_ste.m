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


i2D=3;
mmvet=[0 1];
iLP=0;

mmvet=[3 2]+1;
numodiacc=2;     % numero di modi da accoppiare (sopra, sotto)

mmvet=[3 2];
mmvet=[3];
numodiacc=1;     % numero di modi da accoppiare (sopra, sotto)

%mmvet=[8 9];
%numodiacc=5;     % numero di modi da accoppiare (sopra, sotto)
nomes='dp';


iany=1;  % 1: anisotropia planare, 2: confinata (effetto elettro-ottico )
%iany=0;  % 1: anisotropia planare, 2: confinata (effetto elettro-ottico )
ianys=0;   % strain
iany=0;

if iLP==1
 ianys=0;
 iany=0;
end

ifp=0;
ifp=-5;
%ifp=-10;
ifp=-3;

Frisi0=0.5e-3;
Frisu0=0.9e-3;
Ndisp0=5;
%nk1max=20;
nk1max=30;
alim=0.25;
idyn=0;
iraff=1;
iraff=0;


Frisi0=0.45e-3;
Frisu0=0.65e-3;
Ndisp0=7;
nk1max=40;
alim=0.3;
numodiacc=4;
%%%%%%%%% max vettoriale 250 pt.

%Frisi0=0.45e-3;
%Frisu0=0.95e-3;
%Ndisp0=3;
%nk1max=10;
%alim=0.2;
%numodiacc=2;

ifp=-10;
ifp=-4;

ipolar=2;  % -1, 1, 2 (entrambe)
%ipolar=-1;

fil_str='stein.str';
lambda0=847;   %nm

%% ifp:  -10 alcuni stop campi
%% ifp:  -5  no display campi e no stop
%% ifp:  -4  display campi e no stop

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% calcolo e plot modi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%xroI1=linspace(0,3,30);
%xroI2=linspace(3,6,10);
%xroI=[xroI1 xroI2(2:length(xroI2))];

ploma=2*r_pil;
Np=100;
npfi0=25;
r=linspace(0,ploma,Np);
xroI=r;
xro=r;



global igamveb igamveu GGbext GGuext
global ilo

igamveb=1;   %=1, gamma vero,  =0, riflessione GGext
igamveu=1;
GGbext=1;
GGuext=1;


lfi_inp=2*npfi0;
fimaxi=pi;

lfi_inp=4*npfi0;
fimaxi=2*pi;
iplan=0;

fileeps='Gres_jo2';
Dla=-.005;


disp(' call caloptg '),
nell=3.559124817121504e+000 -6.740211839941768e-004i;

%par_inv=[4.5/2 2.25/2 nell;
%         4.5/2 2.25/2 1;
%         1.8/2 7.2/2 1];
%
%par_inv=[4.5/2  nell;
%         4.5/2  1;
%         1.8/2 1];
%
%par_inv=[];

% variation management (2 nested loops)
% if centers are varying, put a negative number in the corresponding place:
% its absolute value will be the column of the matrix mat_cen where the centers
% must be put.


par_int(:,:,1)=[4.5/2           4.5/2             1.8/2;
                shaf('planar')  shaf('ellipse')   shaf('ellipse');
                2.25/2          2.25/2            7.2/2;
                2.70            2.70              2.70           ;
                shaf('circle')  shaf('circle')    shaf('circle')];

par_int(:,:,2)=[4.5/2           4.5/2             1.8/2;
                shaf('planar')  shaf('ellipse')   shaf('ellipse');
                2.25/2          2.25/2            7.2/2;
                2.70            2.70              2.70           ;
                shaf('circle')  shaf('circle')    shaf('circle')];

par_int(:,:,3)=[4.5/2           4.5/2             1.8/2;
                shaf('planar')  shaf('ellipse')   shaf('ellipse');
                2.25/2          2.25/2            7.2/2;
                2.6             2.6               2.6            ;
                shaf('rhombus') shaf('rhombus')   shaf('rhombus')];

par_int(:,:,4)=[4.5/2           4.5/2             1.8/2;
                 shaf('planar')  shaf('ellipse')   shaf('ellipse');
                 2.25/2          2.25/2            7.2/2;
                 2.6             2.6               2.6            ;
                 shaf('rhombus') shaf('rhombus')   shaf('rhombus')];


ianyv=[0 1 0 1];

icec=1;
cce=[0 5 (2.5+j*5) (-2.5+j*5) -5 (-2.5-j*5) (2.5-j*5)].';
lce=length(cce);
mat_cen(1:lce+1,icec)=[lce; cce];

cce=[0 (2.5+j*5) (-2.5+j*5) (-2.5-j*5) (2.5-j*5)].';
lce=length(cce);
icec=icec+1;
mat_cen(1:lce+1,icec)=[lce; cce];

scas(1,:)='Case 1';
scas(2,:)='Case 2';
scas(3,:)='Case 3';
scas(4,:)='Case 4';
scas(5,:)='Case 5';
scas(6,:)='Case 6';

sp=size(par_int);
spm=sp(end);
if length(sp)==3
 nlooe=sp(end);
 nlooi=sp(end-1);
else
 nlooe=1;
 nlooi=sp(end);
end
nloo=max(nlooi,1);

t00=clock;

for loope=1:nlooe
%for loope=2
 iany=ianyv(loope);
 if spm>0
  if length(sp)==2
   par_inv=par_int(:,loope);
  else
   par_inv=par_int(:,:,loope);
  end
 else
  par_inv=[];
 end
 s=size(par_inv);
 spm=s(2);

 for loop=1:nloo
% for loop=2
   if spm>0
    par_in=par_inv(:,loop);
   else
    par_in=[];
   end
   fiar=find(par_in<0);
   if length(fiar)>0
    cce=mat_cen(:,fix(abs(par_in(fiar))));
    par_in(fiar,:)=cce.';
   end
   par_in

   t0=clock;

   [Eqw,Eout,xro,fian,lambda,delf,gain,ord,nrAzim,Cu,PaDy,ADom,ADcm]=...
   vel(N,T,ro,zeta,xroI,fimaxi,lfi_inp,par_in,r_pil,fil_str,lambda0,...
   ipolar,mmvet,numodiacc,Frisi0,Frisu0,alim,Ndisp0,nk1max,...
   iplan,iraff,idyn,i2D,iLP,iany,ianys,ifp,Dla,fileeps,nomes);



   tim=etime(clock,t0);
   ore=fix(tim/3600);
   minuti=fix((tim-ore*3600)/60);
   secondi=(tim-ore*3600-minuti*60);
   format short
   disp(' Partial elapsed time (hours, minutes, seconds) = ');
   disp([ ore minuti secondi]);

   lamm(:,loop)=lambda./(1+delf);
   demm(:,loop)=delf;
   pamm(:,loop)=ones(size(delf))*loop;
   gamm(:,loop)=gain;
   tymm(:,loop)=PaDy.l';

   sE=size(Eqw);
   Ecmi(1:sE(1),1:sE(2),1:sE(3),loop)=Eqw;
   Eomi.x(1:sE(1),1:sE(2),1:sE(3),loop)=Eout.x;
   Eomi.y(1:sE(1),1:sE(2),1:sE(3),loop)=Eout.y;

   sA=size(ADom);
   Acmi(1:sA(1),1:sA(2),1:sA(3),loop)=ADcm;
   Aomi(1:sA(1),1:sA(2),1:sA(3),loop)=ADom;

   lamme(:,loop,loope)=lambda./(1+delf);
   demme(:,loop,loope)=delf;
   pamme(:,loop,loope)=ones(size(delf))*loop;
   gamme(:,loop,loope)=gain;
   tymme(:,loop,loope)=PaDy.l';
   Ecme(1:sE(1),1:sE(2),1:sE(3),loop,loope)=Eqw;
   Eome.x(1:sE(1),1:sE(2),1:sE(3),loop,loope)=Eout.x;
   Eome.y(1:sE(1),1:sE(2),1:sE(3),loop,loope)=Eout.y;
   Acme(1:sA(1),1:sA(2),1:sA(3),loop,loope)=ADcm;
   Aome(1:sA(1),1:sA(2),1:sA(3),loop,loope)=ADom;

   save longi lamm demm pamm gamm tymm Ecmi Eomi Acmi Aomi
   save longe lamme demme pamme gamme tymme Ecme Eome Acme Aome
%   keyboard
 end

 fiX=find(tymm>0);
 fiY=find(tymm<0);
 coGH=3e5/lambda;
 demg=demm-min(min(demm));
 fic=figure;
 icfig(loope,1)=fic;
 plot(demg(fiX)*coGH,gamm(fiX),'rx'), hold on,
 plot(demg(fiY)*coGH,gamm(fiY),'go'), grid,
 title(['red x: xpol       green o: ypol      ' scas(loope,:)])
 xlabel(' Res. freq  (GHz)')
 ylabel(' Gain  (1/s)')
 Bir=(demm(fiY)-demm(fiX))*coGH;
 Dic=(gamm(fiY)-gamm(fiX));
 tab(:,1)=pamm(fiX);
 tab(:,2)=Bir;
 tab(:,3)=Dic;
 tab
 bir=tab(:,2);
 dic=tab(:,3);

 fic=figure;
 icfig(loope,2)=fic;
 plot(bir,dic,'ro'), grid
 xlabel(' fy-fx  (GHz)')
 ylabel(' Gy-Gx  (1/s)')
 title(scas(loope,:))
 famo=(bir/min(bir)-1)/10;
 famol=1+max(famo);
  for im=1:length(bir);
   gmod=dic(im);
   lmod=bir(im);
   sym=[num2str(im)];
   ht=text(lmod*famol,gmod,sym); set(ht,'fontsize',8);
  end %im

end

   tim=etime(clock,t00);
   ore=fix(tim/3600);
   minuti=fix((tim-ore*3600)/60);
   secondi=(tim-ore*3600-minuti*60);
   format short
   disp(' Total elapsed time (hours, minutes, seconds) = ');
   disp([ ore minuti secondi]);

for kf=1:length(icfig)
 h=figure(icfig(kf,1));
 sh=320*(kf-1);
 pograp=[20+sh 550 300 350];
 set(h,'Position',pograp);
end
for kf=1:length(icfig)
 h=figure(icfig(kf,2));
 sh=320*(kf-1);
 pograp=[20+sh 120 300 350];
 set(h,'Position',pograp);
end
nf=figure;
close(nf)
for kf=1:nf-1
 fig=find(icfig==kf);
 if length(fig)==0
  figure(kf)
 end
end

disp(' Bir in '), icfig
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
   sym=[sum2str(inua) '-' fchar(inur)];
   ht=text(lmod*famol,gmod,sym); set(ht,'fontsize',8);
end %im

imi=min(ord(1,:));
ima=max(ord(1,:));
j1=figure;
j2=figure;
D=((sum(Eqw)));
Dr=reshape(D,length(fian),length(gain));
for im=imi:ima
 fi=find(ord(1,:)==im);
 lf=length(fi);
 if i2D==3
  [du,ip1]=max(Dr(:,fi));
  [du,ip2]=min(Dr(:,fi));
  Ein1=[];
  Ein2=[];
  Eou1=[];
  Eou2=[];
  for k=1:length(ip1)
   Ein1=[Ein1 abs(reshape(Eqw(:,ip1(k),fi(k)),length(xro),1))];
   Ein2=[Ein2 abs(reshape(Eqw(:,ip2(k),fi(k)),length(xro),1))];
%   Eou1=[Eou1 abs(reshape(Eout(:,ip1(k),fi(k)),length(xro),1))];
%   Eou2=[Eou2 abs(reshape(Eout(:,ip2(k),fi(k)),length(xro),1))];
  end
  figure(j1), plot(xro,Ein1,co(:,im+1),xro,Ein2,co1(:,im+1)),
  hold on
%  figure(j2), plot(xro,Eou1,co(:,im+1),xro,Eou2,co1(:,im+1)),
%  hold on
 else
  figure(j1), plot(xro,abs(Eqw(:,fi)),co(:,im+1)), hold on
%  figure(j2), plot(xro,abs(Eout(:,fi)),co(:,im+1)), hold on
 end
end

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
     stri=[' gain = ',num2str(gain(k),7)];
     text(-aax,1.4*aax,stri);
     stri=[' freq = ',num2str(delf(k),7)];
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