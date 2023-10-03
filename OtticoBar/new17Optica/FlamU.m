function lam=Flam(L_in,n_i,iat,la0,g0);


L_i=L_in/1000;
n=pi*50;
n=3;

la=la0;
k0=2*pi/la;

%'contr na FlamU', keyboard
GA=0.5*g0*1e-4/k0;

be=j*k0*n;

n1=n_i(1);
n3=n_i(end);

ns=n_i(2:iat-1);
nd=n_i(iat+1:end-1);
Ls=L_i(2:iat-1);
Ld=L_i(iat+1:end-1);

na=n_i(iat)+j*GA;
La=L_i(iat);


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
lam=2*lamf;


%'in FlamU ', keyboard



