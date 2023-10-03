function [lam,g_gamma]=Flam_gamL(L_in,n_i,iat,icav,la0);

n=3;

L_i=L_in/1000;

la=la0;
k0=2*pi/la;


be=j*k0*n;


na=n_i(iat);
La=L_i(iat);


n1=n_i(1);
n3=n_i(end);

ns=n_i(2:iat-1);
nd=n_i(iat+1:end-1);


Ls=L_i(2:iat-1);
Ld=L_i(iat+1:end-1);

Luv=flipud(Ls);
nuv=flipud(ns);
kk=0;
 [Gs]=gaperdm(kk,0,la,[],[],[],[],0,n1,n,0,Luv,nuv,na);
 Luv=Ld;
 nuv=nd;
 [Gd]=gaperdm(kk,0,la,[],[],[],[],0,n3,n,0,Luv,nuv,na);


%Gafil= 2*La/(La+2*.14);
%'prova ga',

kc=k0*ns(end);
Lc=sum(L_i(icav));
%ze=linspace(0,Lc,100)-Lc/2;
%figure, plot(ze,cos(kc*ze).^2), pausak
Iat=La/2+sin(kc*La)/(2*kc);
Ito=Lc/2+sin(kc*Lc)/(2*kc);
Gat=Iat/Ito;
Gat=2*La/Lc;


g0=-Gat*1e4/(2*La)*log(abs(Gs*Gd));
GA=0.5*g0*1e-4/k0;
na=n_i(iat)+j*GA;
g_gamma=g0;
%'contr na gamL', keyboard
%na=n_i(iat);


n2=ns(1);
nu=nd(end);

Ga1=-(n2-n)/(n2+n);
t1=sqrt(1-Ga1^2);
Van1=[1 Ga1 ; Ga1 1]/t1;
Van1i=[1 -Ga1 ; -Ga1 1]/t1;

Ga2=-(nu-n)/(nu+n);
t2=sqrt(1-Ga2^2);
Van2=[1 Ga2 ; Ga2 1]/t2;
Van2i=[1 -Ga2 ; -Ga2 1]/t2;


T=eye(2);
nx=ns;
Lx=Ls;
for in=1:length(nx)
 ni=nx(in);
 Li=Lx(in);
 del=(ni^2-n^2)/(2*n^2);
 M=[-(1+del) -del; del (1+del)];
 Mv=Li*be*M;
 Ti=expm(Mv); 
 T=Ti*T;
end
Ts=T;

T=eye(2);
nx=nd;
Lx=Ld;
for in=1:length(nx)
 ni=nx(in);
 Li=Lx(in);
 del=(ni^2-n^2)/(2*n^2);
 M=[-(1+del) -del; del (1+del)];
 Mv=Li*be*M;
 Ti=expm(Mv); 
 T=Ti*T;
end
Td=T;

 ni=na;
 Li=La;
 del=(ni^2-n^2)/(n^2);
 Ma=[-(del/2) -del/2; del/2 (del/2)]*be*Li;
 M=[-(1+del/2) -del/2; del/2 (1+del/2)];
 Mv=Li*be*M;
 Tan=expm(Mv); 

Tai=Tan*[-1 -1; 1 1]/2;

%Tn=Van2i*Td*Tan*Ts*Van1;
%Ti=Van2i*Td*Tai*Ts*Van1;

Tn=Td*Tan*Ts;
Ti=Td*Tai*Ts;

% per questi gamma, sempre valutarli da interno a esterno
% per questo Gai ha segno n_i scambiati

%Gas=-(n1-n2)/(n2+n1);
%Gad=-(n3-nu)/(n3+nu);

Gas=-(n1-n)/(n1+n);
Gad=-(n3-n)/(n3+n);

Mn=Gad*Tn(1,1)*Gas+Gad*Tn(1,2)-Tn(2,1)*Gas-Tn(2,2);
Mi=Gad*Ti(1,1)*Gas+Gad*Ti(1,2)-Ti(2,1)*Gas-Ti(2,2);
P=inv(Mi)*Mn;
lamf=P*n/(2*na*La*1e-4);
lam=2*lamf+g0;


%'in Flam_gamu ', keyboard



