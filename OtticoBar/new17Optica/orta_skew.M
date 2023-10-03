%  r_out
%       ____    n4   ____  ^
%       | n2|  n1 |     |      | thick
% __ _|     |___ |     |      v
%         d2   d1
%            n3
% r_in

function [Ge,Gm,Te,Tm,Ge1,Gm1]=orta_skew(phig,tetag,r_in,r_out,r1,r2,d1i,d2i,thicki,lambdai,itetm,Nmodi,ifp);

%clear global
global HE teta err numbet c mu0 divis spp spd sps lp ls ld d d1 d2 n1 n2 n3 n4 eps1 eps2
%global teta c mu0 divis spp spd sps lp ls ld d d1 d2 n1 n2 n3 n4 eps1 eps2

segem=-1;
%segem=1;

if nargin<12
 ifp=-4;
end 

Ge=zeros(1,2);
Gm=Ge;
Ge1=Ge;
Gm1=Gm;
Te=Ge;
Tm=Ge;
d1=d1i*1000;
d2=d2i*1000;
lambdavet=lambdai*1000;
thick=thicki*1000;

theta=tetag*pi/180;
phi=phig*pi/180;

nef=0;

Nstratper=length(thick); % numero di strati in ogni dente
Nlambda=length(lambdai); % numero di strati in ogni dente


n1=r1;
n2=r2;
n3 = r_in;
n4 = r_out;

if length(thick)~=length(n2)
 'errore parametri in orta_gen: thick non consistente con n1,n2',  keyboard
end

d=d1+d2;


if nargin<=9
 itetm=3;
end
if nargin<=10
 numbet=21;
else 
 numbet=Nmodi;
end

flagPlot='Plot_NO';
passo=1;
Clight=2.99792458e8; %nm / ns, quindi frequenze in GHz
mu0=4*pi*1e-7;
eps0=1./(mu0*Clight^2);
NmodiTETM=2*Nmodi;



polariz_inc='TE';
suborta_skew
%S12TETEtotale
Ge(1)= S11TETEtotale;
Ge(2)= S11TETMtotale;
Te(1)= S21TETEtotale;
Te(2)= S21TETMtotale;

Gm(1)= S11TMTMtotale;
Gm(2)= S11TMTEtotale;
Tm(1)= S21TMTMtotale;
Tm(2)= S21TMTEtotale;

Gm1(1)= S22TMTMtotale;
Gm1(2)= S22TMTEtotale;
Ge1(1)= S22TETEtotale;
Ge1(2)= S22TETMtotale;

