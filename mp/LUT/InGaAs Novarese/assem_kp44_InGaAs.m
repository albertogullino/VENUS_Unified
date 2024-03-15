% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% =========================================================================
% This function assembles the 8x8 zinc-blende kdotp Hamiltonian according t
% o Veprek_2007 (which is consistent to Foreman_1997)
% =========================================================================
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Filosofia dietro questa funzione:
%
% 1) Si assembla la matrice Veprek 8x8, dal momento che essa permette di
% implementare operator ordering. Questo è fondamentale perché fa sparire
% soluzioni assurde, addirittura livelli delle sottobande in mezzo al gap.
% 
% 2) Si effettua un cambio di base tra Veprek 8x8 e Chuang 8x8. Il cambio
% di base è riportato tra gli script di Francesco, Veprek_to_Chuang ,
% analitico. Si trova tra gli script della tesi di Cimbri 2017. Notare che
% questa va applicata alle singole matrici di Veprek e non solo alla somma:
% sia a H0, sia a H1L, H1R, H2. Questo perché H1L+H1R ha elementi che vanno
% a 0, ma le singole H1L e H1R no, e questo è fondamentale per l'operator
% ordering.
%
% 3) La formulazione Chuang è necessaria perché permette di essere
% approssimata mediante approssimazione assiale (non è chiaro come la cosa
% sarebbe possibile con la formulazione di Veprek) e poi, con un secondo
% cambio di base (che si trova qua programmato) passare al 4x4.

function [H0,H1L,H1R,H2,Egtot,mesh]=assem_kp44_InGaAs(mesh,kx,ky,xmol,strain,mode)
%
% Conversion to 1/m from angstrom
kx=kx*1e10;
ky=ky*1e10;
%
% Elementary constants
qel=1.6021766208e-019; % Elementary charge, C
h=6.626070040e-34; % Planck constant, J*s
hbar=h/(2*pi); % Reduced Planck constant, J*s
m0=9.10938188E-31; % Electron mass, kg
h2m0=hbar^2/(2*m0);
SR2 = sqrt(2); SR3 = sqrt(3); SR32 = sqrt(3/2);
ind1=[1:4];
ind2=[5:8];
%
% =========================================================================
% GaAs hamiltonian
% =========================================================================
% Material parameters
 Eg = mesh.EgGaAs*qel;
%  Eg = 1.412*qel; % Energy gap, fitted from Gerlach PL data, J
%  Eg= 1.311647647647648*qel;
% %Eg=1.337089089089089*qel; %980nm
% Eg=1.297513513513513*qel;% %0.31 xmol
 Eg_GaAs = 1.412; % Energy gap, fitted from Gerlach PL data, J
%Eg=0;

mc = 0.067; % Electron effective mass, 1995Chuang, p. 709
Delta = 0.34*qel; % spin-orbit coupling, 1995Chuang, p. 709
Ep = 25.7*qel; % Optical matrix parameter, 1995Chuang, p. 709
mc=1/(1/mc-Ep*(Eg+2/3*Delta)/(Eg*(Eg+Delta))); % 2012QiaoPF_OE, eq. (7)
gamma1L = 6.85; % Luttinger parameter, 1995Chuang, p. 709
gamma2L= 2.10; % Luttinger parameter, 1995Chuang, p. 709
gamma3L = 2.90; % Luttinger parameter, 1995Chuang, p. 709
gamma1= gamma1L-Ep/(3*Eg+Delta); % 2012QiaoPF_OE, eq. (8)
gamma2 = gamma2L-0.5*Ep/(3*Eg+Delta); % 2012QiaoPF_OE, eq. (8)
gamma3 = gamma3L-0.5*Ep/(3*Eg+Delta); % 2012QiaoPF_OE, eq. (8)
gamma_bar = (gamma2+gamma3)/2;
C12=5.376;  %Zwang dyn/cm^2
C11=11.879; ;  %Zwang dyn/cm^2
ac=-7.17*qel; % eV
av=1.16 *qel; %eV
b=-1.7*qel; %eV
exx=strain;
eyy=exx;
ezz=-2*(C12/C11)*exx;
%
% -------------------------------------------------------------------------
Astr=ac*(exx+eyy+ezz);
Pstr=-av*(exx+eyy+ezz);
Qstr=-b/2*(exx+eyy-2*ezz);
% Astr=0;
% Pstr=0;
% Qstr=0;

st=b*(C11+2*C12)/C11*strain;
 st=0;
Hs(1,:)=[st 0 0 0];
Hs(2,:)=[ 0 -st 0 0];
Hs(3,:)=[0 0 -st 0];
Hs(4,:)=[ 0 0 0 st];
kz = 0;
kt=sqrt(kx^2+ky^2);
A  = h2m0/mc*(kt^2+kz^2)+Astr; 
P0=sqrt(h2m0*Ep);
P = h2m0*gamma1*(kt^2+kz^2)+Pstr;
Q = h2m0*gamma2*(kt^2-2*kz^2)+Qstr;
% R = h2m0*SR3*(-gamma2*(kx^2-ky^2)+2*1i*gamma3*kx*ky);
% axial approximation
R = h2m0*SR3*(-gamma_bar*(kx^2-ky^2)+2*1i*gamma_bar*kx*ky);
S = h2m0*2*SR3*gamma3*(kx-1i*ky)*kz;
V = 1/sqrt(6)*P0*(kx+1i*ky);
U = 1/SR3*P0*kz;
%
H0_GaAs(1,:) = [Eg+A -SR3*V SR2*U U 0 0 conj(V) SR2*conj(V)];
H0_GaAs(2,:) = [-SR3*conj(V) -P-Q S S/SR2 0 0 -R -SR2*R];
H0_GaAs(3,:) = [SR2*U conj(S) -P+Q SR2*Q -conj(V) -R 0 -SR32*S];
H0_GaAs(4,:) = [U conj(S)/SR2 SR2*Q -P-Delta SR2*conj(V) SR2*R -SR32*S 0];
H0_GaAs(5,:) = [0 0 -V SR2*V Eg+A SR3*conj(V) SR2*U -U];
H0_GaAs(6,:) = [0 0 -conj(R) SR2*conj(R) SR3*V -P-Q -conj(S) conj(S)/SR2];
H0_GaAs(7,:) = [V -conj(R) 0 -SR32*conj(S) SR2*U -S -P+Q -SR2*Q];
H0_GaAs(8,:) = [SR2*V -SR2*conj(R) -SR32*conj(S) 0 -U S/SR2 -SR2*Q -P-Delta];
%
% H1L ---------------------------------------------------------------------
kz = 1;
A = 0;
P = 0;
Q = 0;
R = 0;
M = -h2m0*(gamma1-2*gamma2); % 2007Veprek_PRB, eq. (18)
NL = M - h2m0;
SL = - NL/SR3*(kx-1i*ky)*kz;
NR = - h2m0*6*gamma3 - NL; 
SR = - NR/SR3*(kx-1i*ky)*kz;
V = 0;
U = 1/SR3*P0*kz;
%
H1_GaAsL(1,:) =[                            0,                               0,                                          0,                                          0,                            0,                               0,                                              0,                                              0];
H1_GaAsL(2,:) =[                            0,                               0,            -(3^(1/2)*NL*kz*(kx - ky*1i))/3,            -(6^(1/2)*NL*kz*(kx - ky*1i))/6,                            0,                               0,                                              0,                                              0];
H1_GaAsL(3,:) =[            (6^(1/2)*P0*kz)/3, -(3^(1/2)*NR*kz*(kx + ky*1i))/3,                                          0,                                          0,                            0,                               0,               (kz*(NL - NR)*(ky + kx*1i)*1i)/3, -(2^(1/2)*kz*(ky + kx*1i)*(3*NL + 6*NR)*1i)/18];
H1_GaAsL(4,:) =[            (3^(1/2)*P0*kz)/3, -(6^(1/2)*NR*kz*(kx + ky*1i))/6,                                          0,                                          0,                            0,                               0, -(2^(1/2)*kz*(ky + kx*1i)*(6*NL + 3*NR)*1i)/18,               (kz*(NL - NR)*(ky + kx*1i)*1i)/3];
H1_GaAsL(5,:) =[                            0,                               0,                                          0,                                          0,                            0,                               0,                                              0,                                              0];
H1_GaAsL(6,:) =[                            0,                               0,                                          0,                                          0,                            0,                               0,                 (3^(1/2)*NL*kz*(kx + ky*1i))/3,                -(6^(1/2)*NL*kz*(kx + ky*1i))/6];
H1_GaAsL(7,:) =[                            0,                               0,              (kz*(NL - NR)*(kx + ky*1i))/3, (2^(1/2)*kz*(kx + ky*1i)*(3*NL + 6*NR))/18,            (6^(1/2)*P0*kz)/3,  (3^(1/2)*NR*kz*(kx - ky*1i))/3,                                              0,                                              0];
H1_GaAsL(8,:) =[                            0,                               0, (2^(1/2)*kz*(kx + ky*1i)*(6*NL + 3*NR))/18,              (kz*(NL - NR)*(kx + ky*1i))/3,           -(3^(1/2)*P0*kz)/3, -(6^(1/2)*NR*kz*(kx - ky*1i))/6,                                              0,                                              0];
%
% H1R ---------------------------------------------------------------------
H1_GaAsR(1,:) =[                            0,                               0,                          (6^(1/2)*P0*kz)/3,                          (3^(1/2)*P0*kz)/3,                            0,                               0,                                              0,                                              0];
H1_GaAsR(2,:) =[                            0,                               0,            -(3^(1/2)*NR*kz*(kx - ky*1i))/3,            -(6^(1/2)*NR*kz*(kx - ky*1i))/6,                            0,                               0,                                              0,                                              0];
H1_GaAsR(3,:) =[                            0, -(3^(1/2)*NL*kz*(kx + ky*1i))/3,                                          0,                                          0,                            0,                               0,              -(kz*(NL - NR)*(ky + kx*1i)*1i)/3, -(2^(1/2)*kz*(ky + kx*1i)*(6*NL + 3*NR)*1i)/18];
H1_GaAsR(4,:) =[                            0, -(6^(1/2)*NL*kz*(kx + ky*1i))/6,                                          0,                                          0,                            0,                               0, -(2^(1/2)*kz*(ky + kx*1i)*(3*NL + 6*NR)*1i)/18,              -(kz*(NL - NR)*(ky + kx*1i)*1i)/3];
H1_GaAsR(5,:) =[                            0,                               0,                                          0,                                          0,                            0,                               0,                              (6^(1/2)*P0*kz)/3,                             -(3^(1/2)*P0*kz)/3];
H1_GaAsR(6,:) =[                            0,                               0,                                          0,                                          0,                            0,                               0,                 (3^(1/2)*NR*kz*(kx + ky*1i))/3,                -(6^(1/2)*NR*kz*(kx + ky*1i))/6];
H1_GaAsR(7,:) =[                            0,                               0,             -(kz*(NL - NR)*(kx + ky*1i))/3, (2^(1/2)*kz*(kx + ky*1i)*(6*NL + 3*NR))/18,                            0,  (3^(1/2)*NL*kz*(kx - ky*1i))/3,                                              0,                                              0];
H1_GaAsR(8,:) =[                            0,                               0, (2^(1/2)*kz*(kx + ky*1i)*(3*NL + 6*NR))/18,             -(kz*(NL - NR)*(kx + ky*1i))/3,                            0, -(6^(1/2)*NL*kz*(kx - ky*1i))/6,                                              0,                                              0];
%
% H2 ----------------------------------------------------------------------
kz = 1;
A = h2m0/mc * (kz^2); 
P = h2m0*gamma1*(kz^2);
Q = h2m0*gamma2*(-2*kz^2);
R = 0;
S = 0;
V = 0;
U = 0;
%
H2_GaAs(1,:) = [A -SR3*V SR2*U U 0 0 conj(V) SR2*conj(V)];
H2_GaAs(2,:) = [-SR3*conj(V) -P-Q S S/SR2 0 0 -R -SR2*R];
H2_GaAs(3,:) = [SR2*U conj(S) -P+Q SR2*Q -conj(V) -R 0 -SR32*S];
H2_GaAs(4,:) = [U conj(S)/SR2 SR2*Q -P SR2*conj(V) SR2*R -SR32*S 0];
H2_GaAs(5,:) = [0 0 -V SR2*V A SR3*conj(V) SR2*U -U];
H2_GaAs(6,:) = [0 0 -conj(R) SR2*conj(R) SR3*V -P-Q -conj(S) conj(S)/SR2];
H2_GaAs(7,:) = [V -conj(R) 0 -SR32*conj(S) SR2*U -S -P+Q -SR2*Q];
H2_GaAs(8,:) = [SR2*V -SR2*conj(R) -SR32*conj(S) 0 -U S/SR2 -SR2*Q -P];
% -------------------------------------------------------------------------
H0_GaAs  = H0_GaAs/qel;
H1_GaAsL = H1_GaAsL/qel;
H1_GaAsR = H1_GaAsR/qel;
H2_GaAs  = H2_GaAs/qel;
%
U = zeros(8,8);
if(ky==0), phi = 0; else phi = atan(ky/kx); end
U(1,1) =     exp(-1i*1/2*phi); U(1,5) =   1i*exp(1i*1/2*phi);
U(2,2) =     exp(-1i*3/2*phi); U(2,6) = - 1i*exp(1i*3/2*phi);
U(3,3) =  1i*exp(-1i*1/2*phi); U(3,7) = -    exp(1i*1/2*phi);
U(4,4) = -1i*exp(-1i*1/2*phi); U(4,8) = -    exp(1i*1/2*phi);
U(5,1) =     exp(-1i*1/2*phi); U(5,5) = - 1i*exp(1i*1/2*phi);
U(6,2) =     exp(-1i*3/2*phi); U(6,6) =   1i*exp(1i*3/2*phi);
U(7,3) = -1i*exp(-1i*1/2*phi); U(7,7) = -    exp(1i*1/2*phi);
U(8,4) =  1i*exp(-1i*1/2*phi); U(8,8) = -    exp(1i*1/2*phi);
U = U *1/SR2; 
%
H0_GaAs =  (inv(U.')*(H0_GaAs) *(U.')); 
H2_GaAs =  (inv(U.')*(H2_GaAs) *(U.'));
H1_GaAsL = (inv(U.')*(H1_GaAsL)*(U.')); 
H1_GaAsR = (inv(U.')*(H1_GaAsR)*(U.'));

if mode.fgain==0
H0_GaAs  = H0_GaAs(ind1,ind1)+Hs/qel; 
H2_GaAs  = H2_GaAs(ind1,ind1); 
H1_GaAsL = H1_GaAsL(ind1,ind1); 
 H1_GaAsR = H1_GaAsR(ind1,ind1); 
else
H0_GaAs  = H0_GaAs(ind2,ind2)+Hs/qel; 
H2_GaAs  = H2_GaAs(ind2,ind2); 
H1_GaAsL = H1_GaAsL(ind2,ind2); 
 H1_GaAsR = H1_GaAsR(ind2,ind2); 
end
%
% =========================================================================
% InAs hamiltonian
% =========================================================================
% Material parameters
 Eg = mesh.EgInAs*qel;
%  Eg = 0.354*qel; % Energy gap of InAs chuang
%  Eg=0.328840840840841*qel;
%  % Eg=0.335219219219219*qel; 980nm
%  Eg=0.325297297297297*qel; %0.31 xmol
Eg_InAs = 0.354; % Energy gap of InAs chuang
%Eg=0;
mc = 0.0224; % Electron effective mass, 1995Chuang, p. 709
Delta = 0.38*qel; % spin-orbit coupling, 1995Chuang, p. 709
Ep = 22.2*qel; % Optical matrix parameter, 1995Chuang, p. 709
gamma1L=16.5;% Kim 2009
gamma2L=6.5;% Kim 2009
gamma3L=7.2;% Kim 2009
% modified to take into account 8x8 H chuang 2012
gamma1=gamma1L-Ep/(3*Eg+Delta);% Chuang 2012
gamma2 =gamma2L-0.5*Ep/(3*Eg+Delta); % % Chuang 2012
gamma3 = gamma3L-0.5*Ep/(3*Eg+Delta); % % Chuang 2012
% se fatto due volte per qualche strano motivo viene qualcosa di sensato,
% ma parametri di lutt diventano negativi
mc=1/(1/mc-Ep*(Eg+2/3*Delta)/(Eg*(Eg+Delta))); % 2012QiaoPF_OE, eq. (7)
gamma_bar = (gamma2+gamma3)/2;
C12=4.526;  %Zwang dyn/cm^2
C11=12.329; ;  %Zwang dyn/cm^2
ac=-5.08*qel; % eV
av=2*qel; %eV
b=-1.8*qel ;%eV
exx=strain;
eyy=exx;
ezz=-2*(C12/C11)*exx;
%
% -------------------------------------------------------------------------
Astr=ac*(exx+eyy+ezz);
Pstr=-av*(exx+eyy+ezz);
Qstr=-b/2*(exx+eyy-2*ezz);
% Astr=0;
% Pstr=0;
% Qstr=0;

st=b*(C11+2*C12)/C11*strain;
st=0;
Hs(1,:)=[st 0 0 0];
Hs(2,:)=[ 0 -st 0 0];
Hs(3,:)=[0 0 -st 0];
Hs(4,:)=[ 0 0 0 st];
kz = 0;
kt=sqrt(kx^2+ky^2);
A  = h2m0/mc*(kt^2+kz^2)+Astr; 
P0=sqrt(h2m0*Ep);
P = h2m0*gamma1*(kt^2+kz^2)+Pstr;
Q = h2m0*gamma2*(kt^2-2*kz^2)+Qstr;
% R = h2m0*SR3*(-gamma2*(kx^2-ky^2)+2*1i*gamma3*kx*ky);
% axial approximation
R = h2m0*SR3*(-gamma_bar*(kx^2-ky^2)+2*1i*gamma_bar*kx*ky);
S = h2m0*2*SR3*gamma3*(kx-1i*ky)*kz;
V = 1/sqrt(6)*P0*(kx+1i*ky);
U = 1/SR3*P0*kz;
%
H0_InAs(1,:) = [Eg+A -SR3*V SR2*U U 0 0 conj(V) SR2*conj(V)];
H0_InAs(2,:) = [-SR3*conj(V) -P-Q S S/SR2 0 0 -R -SR2*R];
H0_InAs(3,:) = [SR2*U conj(S) -P+Q SR2*Q -conj(V) -R 0 -SR32*S];
H0_InAs(4,:) = [U conj(S)/SR2 SR2*Q -P-Delta SR2*conj(V) SR2*R -SR32*S 0];
H0_InAs(5,:) = [0 0 -V SR2*V Eg+A SR3*conj(V) SR2*U -U];
H0_InAs(6,:) = [0 0 -conj(R) SR2*conj(R) SR3*V -P-Q -conj(S) conj(S)/SR2];
H0_InAs(7,:) = [V -conj(R) 0 -SR32*conj(S) SR2*U -S -P+Q -SR2*Q];
H0_InAs(8,:) = [SR2*V -SR2*conj(R) -SR32*conj(S) 0 -U S/SR2 -SR2*Q -P-Delta];
%
% H1L ---------------------------------------------------------------------
kz = 1;
A = 0;
P = 0;
Q = 0;
R = 0;
M = -h2m0*(gamma1-2*gamma2); % Renormalized Kane's parameter (Veprek_2007)                      
NL = M - h2m0;
SL = - NL/SR3*(kx-1i*ky)*kz;
NR = - h2m0*6*gamma3 - NL; 
SR = - NR/SR3*(kx-1i*ky)*kz;
V = 0;
U = 1/SR3*P0*kz;
%
H1_InAsL(1,:) =[                            0,                               0,                                          0,                                          0,                            0,                               0,                                              0,                                              0];
H1_InAsL(2,:) =[                            0,                               0,            -(3^(1/2)*NL*kz*(kx - ky*1i))/3,            -(6^(1/2)*NL*kz*(kx - ky*1i))/6,                            0,                               0,                                              0,                                              0];
H1_InAsL(3,:) =[            (6^(1/2)*P0*kz)/3, -(3^(1/2)*NR*kz*(kx + ky*1i))/3,                                          0,                                          0,                            0,                               0,               (kz*(NL - NR)*(ky + kx*1i)*1i)/3, -(2^(1/2)*kz*(ky + kx*1i)*(3*NL + 6*NR)*1i)/18];
H1_InAsL(4,:) =[            (3^(1/2)*P0*kz)/3, -(6^(1/2)*NR*kz*(kx + ky*1i))/6,                                          0,                                          0,                            0,                               0, -(2^(1/2)*kz*(ky + kx*1i)*(6*NL + 3*NR)*1i)/18,               (kz*(NL - NR)*(ky + kx*1i)*1i)/3];
H1_InAsL(5,:) =[                            0,                               0,                                          0,                                          0,                            0,                               0,                                              0,                                              0];
H1_InAsL(6,:) =[                            0,                               0,                                          0,                                          0,                            0,                               0,                 (3^(1/2)*NL*kz*(kx + ky*1i))/3,                -(6^(1/2)*NL*kz*(kx + ky*1i))/6];
H1_InAsL(7,:) =[                            0,                               0,              (kz*(NL - NR)*(kx + ky*1i))/3, (2^(1/2)*kz*(kx + ky*1i)*(3*NL + 6*NR))/18,            (6^(1/2)*P0*kz)/3,  (3^(1/2)*NR*kz*(kx - ky*1i))/3,                                              0,                                              0];
H1_InAsL(8,:) =[                            0,                               0, (2^(1/2)*kz*(kx + ky*1i)*(6*NL + 3*NR))/18,              (kz*(NL - NR)*(kx + ky*1i))/3,           -(3^(1/2)*P0*kz)/3, -(6^(1/2)*NR*kz*(kx - ky*1i))/6,                                              0,                                              0];
%
% H1R ---------------------------------------------------------------------
H1_InAsR(1,:) =[                            0,                               0,                          (6^(1/2)*P0*kz)/3,                          (3^(1/2)*P0*kz)/3,                            0,                               0,                                              0,                                              0];
H1_InAsR(2,:) =[                            0,                               0,            -(3^(1/2)*NR*kz*(kx - ky*1i))/3,            -(6^(1/2)*NR*kz*(kx - ky*1i))/6,                            0,                               0,                                              0,                                              0];
H1_InAsR(3,:) =[                            0, -(3^(1/2)*NL*kz*(kx + ky*1i))/3,                                          0,                                          0,                            0,                               0,              -(kz*(NL - NR)*(ky + kx*1i)*1i)/3, -(2^(1/2)*kz*(ky + kx*1i)*(6*NL + 3*NR)*1i)/18];
H1_InAsR(4,:) =[                            0, -(6^(1/2)*NL*kz*(kx + ky*1i))/6,                                          0,                                          0,                            0,                               0, -(2^(1/2)*kz*(ky + kx*1i)*(3*NL + 6*NR)*1i)/18,              -(kz*(NL - NR)*(ky + kx*1i)*1i)/3];
H1_InAsR(5,:) =[                            0,                               0,                                          0,                                          0,                            0,                               0,                              (6^(1/2)*P0*kz)/3,                             -(3^(1/2)*P0*kz)/3];
H1_InAsR(6,:) =[                            0,                               0,                                          0,                                          0,                            0,                               0,                 (3^(1/2)*NR*kz*(kx + ky*1i))/3,                -(6^(1/2)*NR*kz*(kx + ky*1i))/6];
H1_InAsR(7,:) =[                            0,                               0,             -(kz*(NL - NR)*(kx + ky*1i))/3, (2^(1/2)*kz*(kx + ky*1i)*(6*NL + 3*NR))/18,                            0,  (3^(1/2)*NL*kz*(kx - ky*1i))/3,                                              0,                                              0];
H1_InAsR(8,:) =[                            0,                               0, (2^(1/2)*kz*(kx + ky*1i)*(3*NL + 6*NR))/18,             -(kz*(NL - NR)*(kx + ky*1i))/3,                            0, -(6^(1/2)*NL*kz*(kx - ky*1i))/6,                                              0,                                              0];
%
% H2 ----------------------------------------------------------------------
kz = 1;
A = h2m0/mc * (kz^2); 
P = h2m0*gamma1*(kz^2);
Q = h2m0*gamma2*(-2*kz^2);
R = 0;
S = 0;
V = 0;
U = 0;
%
H2_InAs(1,:) = [A -SR3*V SR2*U U 0 0 conj(V) SR2*conj(V)];
H2_InAs(2,:) = [-SR3*conj(V) -P-Q S S/SR2 0 0 -R -SR2*R];
H2_InAs(3,:) = [SR2*U conj(S) -P+Q SR2*Q -conj(V) -R 0 -SR32*S];
H2_InAs(4,:) = [U conj(S)/SR2 SR2*Q -P SR2*conj(V) SR2*R -SR32*S 0];
H2_InAs(5,:) = [0 0 -V SR2*V A SR3*conj(V) SR2*U -U];
H2_InAs(6,:) = [0 0 -conj(R) SR2*conj(R) SR3*V -P-Q -conj(S) conj(S)/SR2];
H2_InAs(7,:) = [V -conj(R) 0 -SR32*conj(S) SR2*U -S -P+Q -SR2*Q];
H2_InAs(8,:) = [SR2*V -SR2*conj(R) -SR32*conj(S) 0 -U S/SR2 -SR2*Q -P];
% -------------------------------------------------------------------------
H0_InAs  = H0_InAs/qel;
H1_InAsL = H1_InAsL/qel;
H1_InAsR = H1_InAsR/qel;
H2_InAs  = H2_InAs/qel;
%
U = zeros(8,8);
if(ky==0), phi = 0; else phi = atan(ky/kx); end
U(1,1) =     exp(-1i*1/2*phi); U(1,5) =   1i*exp(1i*1/2*phi);
U(2,2) =     exp(-1i*3/2*phi); U(2,6) = - 1i*exp(1i*3/2*phi);
U(3,3) =  1i*exp(-1i*1/2*phi); U(3,7) = -    exp(1i*1/2*phi);
U(4,4) = -1i*exp(-1i*1/2*phi); U(4,8) = -    exp(1i*1/2*phi);
U(5,1) =     exp(-1i*1/2*phi); U(5,5) = - 1i*exp(1i*1/2*phi);
U(6,2) =     exp(-1i*3/2*phi); U(6,6) =   1i*exp(1i*3/2*phi);
U(7,3) = -1i*exp(-1i*1/2*phi); U(7,7) = -    exp(1i*1/2*phi);
U(8,4) =  1i*exp(-1i*1/2*phi); U(8,8) = -    exp(1i*1/2*phi);
U = U *1/SR2; 
%
H0_InAs =  (inv(U.')*(H0_InAs) *(U.'));
H2_InAs =  (inv(U.')*(H2_InAs) *(U.')); 
H1_InAsL = (inv(U.')*(H1_InAsL)*(U.')); 
H1_InAsR = (inv(U.')*(H1_InAsR)*(U.')); 
if mode.fgain==0
H0_InAs  = H0_InAs(ind1,ind1)+Hs/qel; 
H2_InAs  = H2_InAs(ind1,ind1); 
H1_InAsL = H1_InAsL(ind1,ind1); 
H1_InAsR = H1_InAsR(ind1,ind1); 
else
H0_InAs  = H0_InAs(ind2,ind2)+Hs/qel; 
H2_InAs  = H2_InAs(ind2,ind2); 
H1_InAsL = H1_InAsL(ind2,ind2); 
H1_InAsR = H1_InAsR(ind2,ind2);
end
%
% =========================================================================
% InGaAs hamiltonian interpolation
% =========================================================================
H0 =  xmol*H0_GaAs  +(1- xmol)*H0_InAs;
H1L =xmol*H1_GaAsL + (1-xmol)*H1_InAsL;
H1R =xmol*H1_GaAsR + (1-xmol)*H1_InAsR;
H2 =  xmol*H2_GaAs  +(1- xmol)*H2_InAs;
Egtot =  xmol*Eg_GaAs  +(1- xmol)*Eg_InAs;
% Eg_InGaAs=0.36 +0.629*xmol + 0.426*xmol^2;
% Egtot=Eg_InGaAs;
% H0(1,1)=H0(1,1)+Eg_InGaAs;
mesh.c1=H0(1,1);
mesh.v1=H0(2,2);
mesh.v2=H0(3,3);


% in H0 termini diagonale e a zero senza strain, aggiungendo strain nuovi termini(plot con leggero e pesante)
