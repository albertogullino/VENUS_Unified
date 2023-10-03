%  r_out
%       ____    n4   ____  ^
%       | n2|  n1 |     |      | thick
% __ _|     |___ |     |      v
%         d2   d1
%            n3
% r_in

function [Ttote,Ttotm]=orta_skewML(phig,tetag,r_in,r_out,d_ini,d_outi,r1,r2,d1i,d2i,thicki,lambdai,Nmodi,ifp);

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

d_in=d_ini*1000;
d_out=d_outi*1000;

lambdavet=lambdai*1000;
thick=thicki*1000;

theta_in=tetag*pi/180;
phi_in=phig*pi/180;

nef=0;

Nstratper=length(thick); % numero di strati in ogni dente
Nlambda=length(lambdai); % numero di strati in ogni dente


n1=r1;
n2=r2;
%n3 = r_in;
%n4 = r_out;

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
if ifp==1
 flagPlot='Plot_SI';
end
passo=1;
Clight=2.99792458e8; %nm / ns, quindi frequenze in GHz
mu0=4*pi*1e-7;
eps0=1./(mu0*Clight^2);
NmodiTETM=2*Nmodi;

%'controllo prima skew', keyboard

suborta_skewML

%'controllo dopo skew', keyboard

S11_in=[SIN11(1,[1 Nmodi+1]); SIN11(Nmodi+1,[1 Nmodi+1])];
S12_in=[SIN12(1,[1 Nmodi+1]); SIN12(Nmodi+1,[1 Nmodi+1])];
S21_in=[SIN21(1,[1 Nmodi+1]); SIN21(Nmodi+1,[1 Nmodi+1])];
S22_in=[SIN22(1,[1 Nmodi+1]); SIN22(Nmodi+1,[1 Nmodi+1])];
        sinv=inv(S12_in);
	T11_in=S21_in-S22_in*sinv*S11_in;
	T12_in=S22_in*sinv;
	T21_in=-sinv*S11_in;
	T22_in=sinv;
	T_in=[T11_in T12_in; T21_in T22_in];



S11_out=[SOUT11(1,[1 Nmodi+1]); SOUT11(Nmodi+1,[1 Nmodi+1])];
S12_out=[SOUT12(1,[1 Nmodi+1]); SOUT12(Nmodi+1,[1 Nmodi+1])];
S21_out=[SOUT21(1,[1 Nmodi+1]); SOUT21(Nmodi+1,[1 Nmodi+1])];
S22_out=[SOUT22(1,[1 Nmodi+1]); SOUT22(Nmodi+1,[1 Nmodi+1])];
        soutv=inv(S12_out);
	T11_out=S21_out-S22_out*soutv*S11_out;
	T12_out=S22_out*soutv;
	T21_out=-soutv*S11_out;
	T22_out=soutv;
	T_out=[T11_out T12_out; T21_out T22_out];
	
S11_t=[St_11(1,[1 Nmodi+1],1); St_11(Nmodi+1,[1 Nmodi+1],1)];
S12_t=[St_12(1,[1 Nmodi+1],1); St_12(Nmodi+1,[1 Nmodi+1],1)];
S21_t=[St_21(1,[1 Nmodi+1],1); St_21(Nmodi+1,[1 Nmodi+1],1)];
S22_t=[St_22(1,[1 Nmodi+1],1); St_22(Nmodi+1,[1 Nmodi+1],1)];	
        soutv=inv(S12_t);
	T11_t=S21_t-S22_t*soutv*S11_t;
	T12_t=S22_t*soutv;
	T21_t=-soutv*S11_t;
	T22_t=soutv;
	T_t=[T11_t T12_t; T21_t T22_t];
	
	T_g=inv(T_out)*T_t*inv(T_in);
	
	


Ttote=[T_g(1,[1 3]); T_g(3,[1 3])];    
Ttotm=[T_g(2,[2 4]); T_g(4,[2 4])];    



%'IN suborta ML 1D', keyboard