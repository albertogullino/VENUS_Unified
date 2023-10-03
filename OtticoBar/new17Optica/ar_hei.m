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
r_pil=8;

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


nomes='dp';


iany=1;  % 1: anisotropia planare, 2: confinata (effetto elettro-ottico )
%iany=0;  % 1: anisotropia planare, 2: confinata (effetto elettro-ottico )
ianys=0;   % strain
iany=0;

iLP=1;
if iLP==1
 ianys=0;
 iany=0;
end

ifp=0;
ifp=-5;
ifp=1;
ifp=-10;
%ifp=-4;

idyn=0;
iraff=0;
nmasce=5;

Dlam_mod(1)=0;   %in nm
Dlam_mod(2)=1;   %in nm

nK_dis=30;
Ev_or_Od='Even';
%Ev_or_Od='Odd';
numodiacc=9;
alim=.1;

numodiacc=5;
nK_dis=20;
%nk1max=50;
%ifp=-4;


ipolar=2;  % -1, 1, 2 (entrambe)
%ipolar=-1;

fil_str='har.str';

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

ploma=1.5*r_pil;
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


% variation management (2 nested loops)
% if centers are varying, put a negative number in the corresponding place:
% its absolute value will be the column of the matrix mat_cen where the centers
% must be put.


%par_ind{1,1}= 3;
%par_ind{2,1}= 3;
%       radii.array{pcdu,2}=ipcen-1;
%       radii.array{pcdu,3}=cemem;
%       radii.array{pcdu,4}=shapes;
%       radii.array{pcdu,5}=Ry;
%       radii.array{pcdu,6}=Rx;
%       radii.array{pcdu,7}=Delta;

par_ind{1,1}= 1;
par_ind{1,3}= [0; 6; 3+5j];
par_ind{1,4}= [4; 2; 1];
par_ind{1,5}= [2; 2; 2.5];
par_ind{1,6}= [2; 2; 2.5];
par_ind{1,7}= [-.2; .4; 0];
par_ind{2,1}= 9;

%par_ind(3,1)= 0;
%par_in(2,1)= shaf('rhombus');


Dv=[1.6];
dperc=[-2 -1];
dox=[9  5];



%ianyv=[0 1 0 1];
npairv=[18 20 22 24 26 27];
ianyv=ones(size(npairv));

icec=1;
fcex=[-3.5 0 3.5];
fcey=[-3.5 0 3.5];
Fcx=ones(size(fcey'))*fcex;
Fcy=fcey'*ones(size(fcex));
Fxy=Fcx+j*Fcy;
cce=reshape(Fxy,prod(size(Fxy)),1);
%figure, plot(real(cce),imag(cce),'o'), pausak

%cce=[0 5 (2.5+j*5) (-2.5+j*5) -5 (-2.5-j*5) (2.5-j*5)].';

lce=length(cce);
mat_cen(1:lce+1,icec)=[lce; cce];
Ar{1}={lce, cce.'};

fcex=[-5.25 -1.75  1.75 5.25];
fcey=[-5.25 -1.75  1.75 5.25];
Fcx=ones(size(fcey'))*fcex;
Fcy=fcey'*ones(size(fcex));
Fxy=Fcx+j*Fcy;
cce=reshape(Fxy,prod(size(Fxy)),1);
%figure, plot(real(cce),imag(cce),'o'), pausak

lce=length(cce);
icec=icec+1;
mat_cen(1:lce+1,icec)=[lce; cce];
Ar{2}={lce, cce.'};

scas(1,:)='Case 1';
scas(2,:)='Case 2';
scas(3,:)='Case 3';
scas(4,:)='Case 4';
scas(5,:)='Case 5';
scas(6,:)='Case 6';

t00=clock;


nloope=length(Dv);
nloopi=length(dperc);
iany=0;

%for loope=1:nloope
for loope=1
 D=Dv(loope);
% par_ind(1,1)=D;

% for loop=1:nloopi
 for loop=1
%  dpe=dperc(loop);
%  par_ind(4,1)=dpe;
%  par_ind(5,1)=dox(loop);
%  for ki=1:length(par_ind)
%   par_in{ki,1}=par_ind(ki,1);
%  end
%   fiar=find(fix(par_ind(:,1))<0);
%   fiar1=find(par_ind(:,1)==5);
%   if length(fiar)*length(fiar1)>0
%    cce=Ar{fix(abs(par_ind(fiar)))};
%    par_in{fiar}=cce;
%   end
  par_in=par_ind
  pausak
%   par_in=[];
   t0=clock;

   fig_ind.a=loope;
   fig_ind.b=loop;

   [Eqw,Eout,xro,fian,lambda,delf,gain,ord,nrAzim,Cu,PaDy,ADom,ADcm]=...
   vel(N,T,ro,zeta,xroI,fimaxi,lfi_inp,par_in,r_pil,fil_str,...
   ipolar,Ev_or_Od,nmasce,numodiacc,Dlam_mod,nK_dis,alim,fig_ind,...
   iplan,iraff,idyn,i2D,iLP,iany,ianys,ifp,Dla,fileeps,nomes);

%   Frisi0,Frisu0,alim,Ndisp0,nk1max,...


   tim=etime(clock,t0);
   ore=fix(tim/3600);
   minuti=fix((tim-ore*3600)/60);
   secondi=(tim-ore*3600-minuti*60);
   format short
   disp(' Partial elapsed time (hours, minutes, seconds) = ');
   disp([ ore minuti secondi]);
 end

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


ick=1;
llpi=length(loopi);
for kl=loopi
 ty=reshape(tymme(:,kl,:),2,nloope);
 fiY=find(ty<0)
 fiX=find(ty>0)
 dG=reshape(gamme(:,kl,:),2,nloope);
 direl(ick,:)=((dG(fiY)-dG(fiX))./mean(dG)')';
 dfr=reshape(demme(:,kl,:),2,nloope);
 bivet(ick,:)=((dfr(fiY)-dfr(fiX))*coGH)';
 ick=ick+1;
end
figure, plot(npairv,direl)
figure, plot(npairv,bivet)
save fihco
