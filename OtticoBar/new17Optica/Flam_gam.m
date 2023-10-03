function [lam,g0]=Flam_gam(L_in,n_i,iat,la0);


L_i=L_in/1000;
n=pi*50;

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
[Gs]=gaperdm(0,0,la,[],[],[],[],0,na,n,0,Luv,nuv,n1);
Luv=Ld;
nuv=nd;
[Gd]=gaperdm(0,0,la,[],[],[],[],0,na,n,0,Luv,nuv,n3);

%Gafil= 2*La/(La+2*.14);
%'prova ga',

kc=k0*ns(end);
Lc=sum(L_i(iat+[-1 0 1]));
%ze=linspace(0,Lc,100)-Lc/2;
%figure, plot(ze,cos(kc*ze).^2), pausak
Iat=La/2+sin(kc*La)/(2*kc);
Ito=Lc/2+sin(kc*Lc)/(2*kc);
Gat=Iat/Ito;
Gat=2*La/Lc;


g0=-Gat*1e4/(2*La)*log(abs(Gs*Gd))
GA=0.5*g0*1e-4/k0;
na=n_i(iat)+j*GA;
%'contr', keyboard


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
 del=(ni^2-n^2)/(2*n^2);
 M=[-(1+del) -del; del (1+del)];
 Mv=Li*be*M;
 Tan=expm(Mv); 

Tai=Tan*[-1 -1; 1 1]/2;

Tn=Van2i*Td*Tan*Ts*Van1;
Ti=Van2i*Td*Tai*Ts*Van1;


% per questi gamma, sempre valutarli da interno a esterno
% per questo Gai ha segno n_i scambiati

Gas=-(n1-n2)/(n2+n1);
%ti=sqrt(1-Gai^2);
%Vai=1/ti*[1 -Gai; -Gai 1];
%Vaii=1/ti*[1 Gai; Gai 1];

Gad=-(n3-nu)/(n3+nu);
%tu=sqrt(1-Gau^2);
%Vau=[1 -Gau ; -Gau 1]/tu;
%Vaui=[1 Gau ; Gau 1]/tu;

Mn=Gad*Tn(1,1)*Gas+Gad*Tn(1,2)-Tn(2,1)*Gas-Tn(2,2);
Mi=Gad*Ti(1,1)*Gas+Gad*Ti(1,2)-Ti(2,1)*Gas-Ti(2,2);
P=inv(Mi)*Mn;
lamf=P*n/(2*na*La*1e-4);
lam=2*lamf+g0;


%'in Flam ', keyboard



