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

function [H0,H1L,H1R,H2]=assem_kp44_AlGaAs(kx,ky,xmol)
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
%
% =========================================================================
% GaAs hamiltonian
% =========================================================================
% Material parameters
Eg = 1.412*qel; % Energy gap, fitted from Gerlach PL data, J
mc = 0.067; % Electron effective mass, 1995Chuang, p. 709
gamma1 = 6.85; % Luttinger parameter, 1995Chuang, p. 709
gamma2 = 2.10; % Luttinger parameter, 1995Chuang, p. 709
gamma3 = 2.90; % Luttinger parameter, 1995Chuang, p. 709
Delta = 0.34*qel; % spin-orbit coupling, 1995Chuang, p. 709
Ep = 25.7*qel; % Optical matrix parameter, 1995Chuang, p. 709
% Modified Luttinger parameters to account for the effect of remote bands
% (only for k.p formulations including both conduction and valence bands)
mc=1/(1/mc-Ep*(Eg+2/3*Delta)/(Eg*(Eg+Delta))); % 2012QiaoPF_OE, eq. (7)
gamma1 = gamma1-Ep/(3*Eg+Delta); % 2012QiaoPF_OE, eq. (8)
gamma2 = gamma2-0.5*Ep/(3*Eg+Delta); % 2012QiaoPF_OE, eq. (8)
gamma3 = gamma3-0.5*Ep/(3*Eg+Delta); % 2012QiaoPF_OE, eq. (8)
gamma_bar = (gamma2+gamma3)/2;
%
% -------------------------------------------------------------------------
kz = 0;
kt=sqrt(kx^2+ky^2);
A  = h2m0/mc*(kt^2+kz^2); 
P0=sqrt(h2m0*Ep);
P = h2m0*gamma1*(kt^2+kz^2);
Q = h2m0*gamma2*(kt^2-2*kz^2);
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
H0_GaAs =  (inv(U.')*(H0_GaAs) *(U.')); H0_GaAs  = H0_GaAs(1:4,1:4); 
H2_GaAs =  (inv(U.')*(H2_GaAs) *(U.')); H2_GaAs  = H2_GaAs(1:4,1:4); 
H1_GaAsL = (inv(U.')*(H1_GaAsL)*(U.')); H1_GaAsL = H1_GaAsL(1:4,1:4); 
H1_GaAsR = (inv(U.')*(H1_GaAsR)*(U.')); H1_GaAsR = H1_GaAsR(1:4,1:4); 
%
% =========================================================================
% AlAs hamiltonian
% =========================================================================
% Material parameters
Eg = 2.6590*qel; % Energy gap of GaAs + DeltaEg (DeltaEg=1.427 eV)
mc = 0.124; % Electron effective mass, 1995Chuang, p. 709
gamma1 = 3.45; % Luttinger parameter, 1995Chuang, p. 709
gamma2 = 0.68; % Luttinger parameter, 1995Chuang, p. 709
gamma3 = 1.29; % Luttinger parameter, 1995Chuang, p. 709
Delta = 0.28*qel; % spin-orbit coupling, 1995Chuang, p. 709
Ep = 21.1*qel; % Optical matrix parameter, 1995Chuang, p. 709
% Modified Luttinger parameters to account for the effect of remote bands
% (only for k.p formulations including both conduction and valence bands)
mc=1/(1/mc-Ep*(Eg+2/3*Delta)/(Eg*(Eg+Delta))); % 2012QiaoPF_OE, eq. (7)
gamma1 = gamma1-Ep/(3*Eg+Delta); % 2012QiaoPF_OE, eq. (8)
gamma2 = gamma2-0.5*Ep/(3*Eg+Delta); % 2012QiaoPF_OE, eq. (8)
gamma3 = gamma3-0.5*Ep/(3*Eg+Delta); % 2012QiaoPF_OE, eq. (8)
gamma_bar = (gamma2+gamma3)/2;
%
% -------------------------------------------------------------------------
kz = 0;
kt=sqrt(kx^2+ky^2);
A  = h2m0/mc*(kt^2+kz^2); 
P0=sqrt(h2m0*Ep);
P = h2m0*gamma1*(kt^2+kz^2);
Q = h2m0*gamma2*(kt^2-2*kz^2);
% R = h2m0*SR3*(-gamma2*(kx^2-ky^2)+2*1i*gamma3*kx*ky);
% axial approximation
R = h2m0*SR3*(-gamma_bar*(kx^2-ky^2)+2*1i*gamma_bar*kx*ky);
S = h2m0*2*SR3*gamma3*(kx-1i*ky)*kz;
V = 1/sqrt(6)*P0*(kx+1i*ky);
U = 1/SR3*P0*kz;
%
H0_AlAs(1,:) = [Eg+A -SR3*V SR2*U U 0 0 conj(V) SR2*conj(V)];
H0_AlAs(2,:) = [-SR3*conj(V) -P-Q S S/SR2 0 0 -R -SR2*R];
H0_AlAs(3,:) = [SR2*U conj(S) -P+Q SR2*Q -conj(V) -R 0 -SR32*S];
H0_AlAs(4,:) = [U conj(S)/SR2 SR2*Q -P-Delta SR2*conj(V) SR2*R -SR32*S 0];
H0_AlAs(5,:) = [0 0 -V SR2*V Eg+A SR3*conj(V) SR2*U -U];
H0_AlAs(6,:) = [0 0 -conj(R) SR2*conj(R) SR3*V -P-Q -conj(S) conj(S)/SR2];
H0_AlAs(7,:) = [V -conj(R) 0 -SR32*conj(S) SR2*U -S -P+Q -SR2*Q];
H0_AlAs(8,:) = [SR2*V -SR2*conj(R) -SR32*conj(S) 0 -U S/SR2 -SR2*Q -P-Delta];
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
H1_AlAsL(1,:) =[                            0,                               0,                                          0,                                          0,                            0,                               0,                                              0,                                              0];
H1_AlAsL(2,:) =[                            0,                               0,            -(3^(1/2)*NL*kz*(kx - ky*1i))/3,            -(6^(1/2)*NL*kz*(kx - ky*1i))/6,                            0,                               0,                                              0,                                              0];
H1_AlAsL(3,:) =[            (6^(1/2)*P0*kz)/3, -(3^(1/2)*NR*kz*(kx + ky*1i))/3,                                          0,                                          0,                            0,                               0,               (kz*(NL - NR)*(ky + kx*1i)*1i)/3, -(2^(1/2)*kz*(ky + kx*1i)*(3*NL + 6*NR)*1i)/18];
H1_AlAsL(4,:) =[            (3^(1/2)*P0*kz)/3, -(6^(1/2)*NR*kz*(kx + ky*1i))/6,                                          0,                                          0,                            0,                               0, -(2^(1/2)*kz*(ky + kx*1i)*(6*NL + 3*NR)*1i)/18,               (kz*(NL - NR)*(ky + kx*1i)*1i)/3];
H1_AlAsL(5,:) =[                            0,                               0,                                          0,                                          0,                            0,                               0,                                              0,                                              0];
H1_AlAsL(6,:) =[                            0,                               0,                                          0,                                          0,                            0,                               0,                 (3^(1/2)*NL*kz*(kx + ky*1i))/3,                -(6^(1/2)*NL*kz*(kx + ky*1i))/6];
H1_AlAsL(7,:) =[                            0,                               0,              (kz*(NL - NR)*(kx + ky*1i))/3, (2^(1/2)*kz*(kx + ky*1i)*(3*NL + 6*NR))/18,            (6^(1/2)*P0*kz)/3,  (3^(1/2)*NR*kz*(kx - ky*1i))/3,                                              0,                                              0];
H1_AlAsL(8,:) =[                            0,                               0, (2^(1/2)*kz*(kx + ky*1i)*(6*NL + 3*NR))/18,              (kz*(NL - NR)*(kx + ky*1i))/3,           -(3^(1/2)*P0*kz)/3, -(6^(1/2)*NR*kz*(kx - ky*1i))/6,                                              0,                                              0];
%
% H1R ---------------------------------------------------------------------
H1_AlAsR(1,:) =[                            0,                               0,                          (6^(1/2)*P0*kz)/3,                          (3^(1/2)*P0*kz)/3,                            0,                               0,                                              0,                                              0];
H1_AlAsR(2,:) =[                            0,                               0,            -(3^(1/2)*NR*kz*(kx - ky*1i))/3,            -(6^(1/2)*NR*kz*(kx - ky*1i))/6,                            0,                               0,                                              0,                                              0];
H1_AlAsR(3,:) =[                            0, -(3^(1/2)*NL*kz*(kx + ky*1i))/3,                                          0,                                          0,                            0,                               0,              -(kz*(NL - NR)*(ky + kx*1i)*1i)/3, -(2^(1/2)*kz*(ky + kx*1i)*(6*NL + 3*NR)*1i)/18];
H1_AlAsR(4,:) =[                            0, -(6^(1/2)*NL*kz*(kx + ky*1i))/6,                                          0,                                          0,                            0,                               0, -(2^(1/2)*kz*(ky + kx*1i)*(3*NL + 6*NR)*1i)/18,              -(kz*(NL - NR)*(ky + kx*1i)*1i)/3];
H1_AlAsR(5,:) =[                            0,                               0,                                          0,                                          0,                            0,                               0,                              (6^(1/2)*P0*kz)/3,                             -(3^(1/2)*P0*kz)/3];
H1_AlAsR(6,:) =[                            0,                               0,                                          0,                                          0,                            0,                               0,                 (3^(1/2)*NR*kz*(kx + ky*1i))/3,                -(6^(1/2)*NR*kz*(kx + ky*1i))/6];
H1_AlAsR(7,:) =[                            0,                               0,             -(kz*(NL - NR)*(kx + ky*1i))/3, (2^(1/2)*kz*(kx + ky*1i)*(6*NL + 3*NR))/18,                            0,  (3^(1/2)*NL*kz*(kx - ky*1i))/3,                                              0,                                              0];
H1_AlAsR(8,:) =[                            0,                               0, (2^(1/2)*kz*(kx + ky*1i)*(3*NL + 6*NR))/18,             -(kz*(NL - NR)*(kx + ky*1i))/3,                            0, -(6^(1/2)*NL*kz*(kx - ky*1i))/6,                                              0,                                              0];
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
H2_AlAs(1,:) = [A -SR3*V SR2*U U 0 0 conj(V) SR2*conj(V)];
H2_AlAs(2,:) = [-SR3*conj(V) -P-Q S S/SR2 0 0 -R -SR2*R];
H2_AlAs(3,:) = [SR2*U conj(S) -P+Q SR2*Q -conj(V) -R 0 -SR32*S];
H2_AlAs(4,:) = [U conj(S)/SR2 SR2*Q -P SR2*conj(V) SR2*R -SR32*S 0];
H2_AlAs(5,:) = [0 0 -V SR2*V A SR3*conj(V) SR2*U -U];
H2_AlAs(6,:) = [0 0 -conj(R) SR2*conj(R) SR3*V -P-Q -conj(S) conj(S)/SR2];
H2_AlAs(7,:) = [V -conj(R) 0 -SR32*conj(S) SR2*U -S -P+Q -SR2*Q];
H2_AlAs(8,:) = [SR2*V -SR2*conj(R) -SR32*conj(S) 0 -U S/SR2 -SR2*Q -P];
% -------------------------------------------------------------------------
H0_AlAs  = H0_AlAs/qel;
H1_AlAsL = H1_AlAsL/qel;
H1_AlAsR = H1_AlAsR/qel;
H2_AlAs  = H2_AlAs/qel;
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
H0_AlAs =  (inv(U.')*(H0_AlAs) *(U.')); H0_AlAs  = H0_AlAs(1:4,1:4); 
H2_AlAs =  (inv(U.')*(H2_AlAs) *(U.')); H2_AlAs  = H2_AlAs(1:4,1:4); 
H1_AlAsL = (inv(U.')*(H1_AlAsL)*(U.')); H1_AlAsL = H1_AlAsL(1:4,1:4); 
H1_AlAsR = (inv(U.')*(H1_AlAsR)*(U.')); H1_AlAsR = H1_AlAsR(1:4,1:4); 
%
% =========================================================================
% AlGaAs hamiltonian interpolation
% =========================================================================
H0 =  (1-xmol)*H0_GaAs  + xmol*H0_AlAs;
H1L = (1-xmol)*H1_GaAsL + xmol*H1_AlAsL;
H1R = (1-xmol)*H1_GaAsR + xmol*H1_AlAsR;
H2 =  (1-xmol)*H2_GaAs  + xmol*H2_AlAs;
