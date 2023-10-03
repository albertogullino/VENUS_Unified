%  r_out
%       ____    n4   ____  ^
%       | n2|  n1 |     |      | thick
% __ _|     |___ |     |      v
%         d2   d1
%            n3
% r_in

function [Ga,Tr]=orta_skewTr(kxi,kyi,r_in,r_out,r1,r2,d1i,d2i,thicki,lambdai,Nmodi,ifp);

%clear global
global HE teta err numbet c mu0 divis spp spd sps lp ls ld d d1 d2 n1 n2 n3 n4 eps1 eps2
%global teta c mu0 divis spp spd sps lp ls ld d d1 d2 n1 n2 n3 n4 eps1 eps2

if nargin<12
 ifp=-4;
end 

Ge=zeros(1,2);
Gm=Ge;
Te=Ge;
Tm=Ge;
d1=d1i*1000;
d2=d2i*1000;
lambdavet=lambdai*1000;
thick=thicki*1000;
kx=kxi/1000;
ky=kyi/1000;


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
%flagPlot='Plot_SI';
passo=1;
Clight=2.99792458e8; %nm / ns, quindi frequenze in GHz
mu0=4*pi*1e-7;
eps0=1./(mu0*Clight^2);
NmodiTETM=2*Nmodi;


%'den', keyboard

polariz_inc='TE';
suborta_skew1
%S12TETEtotale
Ga(1,1)= S11TETEtotale;
Ga(2,1)= S11TMTEtotale;
Ga(2,2)= S11TMTMtotale;
Ga(1,2)= S11TETMtotale;
Tr(1,1)= S21TETEtotale;
Tr(2,1)= S21TMTEtotale;
Tr(2,2)= S21TMTMtotale;
Tr(1,2)= S21TETMtotale;
%'dopo Ga kx ky', keyboard

