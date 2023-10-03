%  r_out
%       ____    n4   ____  ^
%       | n2|  n1 |     |      | thick
% __ _|     |___ |     |      v
%         d2   d1
%            n3
% r_in

function [T11,T22,T12,T21,s11,s12,s21,s22]=orta_skewTOTr(phig,tetag,r_in,r_out,r1,r2,d1i,d2i,thicki,lambdai,Nmodi,ifp,iplotCa,segem);

%clear global
global HE teta err numbet c mu0 divis spp spd sps lp ls ld d d1 d2 n1 n2 n3 n4 eps1 eps2
%global teta c mu0 divis spp spd sps lp ls ld d d1 d2 n1 n2 n3 n4 eps1 eps2

if nargin<12
 ifp=-4;
end 

z11=zeros(2,2);
z12=z11;
z21=z11;
z22=z11;
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
if iplotCa==1
 flagPlot='Plot_SI';
end 
passo=1;

%'in ototr', keyboard
Clight=2.99792458e8; %nm / ns, quindi frequenze in GHz
mu0=4*pi*1e-7;
eps0=1./(mu0*Clight^2);
NmodiTETM=2*Nmodi;


polariz_inc='TE';
polariz_inc='TM';
suborta_skew
%S12TETEtotale
z11(1,1)= S11TETEtotale;
z11(1,2)= S11TETMtotale;
z11(2,1)= S11TMTEtotale;
z11(2,2)= S11TMTMtotale;

z12(1,1)= S12TETEtotale;
z12(1,2)= S12TETMtotale;
z12(2,1)= S12TMTEtotale;
z12(2,2)= S12TMTMtotale;

z21(1,1)= S21TETEtotale;
z21(1,2)= S21TETMtotale;
z21(2,1)= S21TMTEtotale;
z21(2,2)= S21TMTMtotale;

z22(1,1)= S22TETEtotale;
z22(1,2)= S22TETMtotale;
z22(2,1)= S22TMTEtotale;
z22(2,2)= S22TMTMtotale;

%zz=[z11 z12; z21 z22];
        s11=z11;
        s21=z21;
        s12=z12;
        s22=z22;

%        T11=z11;
%        T21=z21;
%        T12=z12;
%        T22=z22;
                
        sinv=inv(s12);
	T11=s21-s22*sinv*s11;
	T12=s22*sinv;
	T21=-sinv*s11;
	T22=sinv;




S =[s11 s12 ; s21 s22];
U=S*S';
%'cont T', keyboard