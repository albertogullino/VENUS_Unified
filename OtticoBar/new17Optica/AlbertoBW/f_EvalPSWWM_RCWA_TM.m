function [kzTM,kTG2,VTM,ITM]=f_EvalPSWWM_RCWA_TM(k0,NModes,csi,ky,E,F)
%
% [kzTM,kTG2,VTM,ITM]=f_EvalPSWWM_RCWA_TM(k0,NModes,csi,ky,E,F)
%
%------------------------------------------------------------------------
% Synthesis of PSWW TM modes
%------------------------------------------------------------------------
% This function calculates the propagation constants (along z) and the
% change of basis matrices from the generic Floquet modes basis to the PSWW
% modes, as the solution of an eigenvalue problem. This eigenvalue problem
% has been solved by using the Lifeng Li formulation, based on the Fourier
% transform inverse rule (in place of the Laurent rule).
%
% NModes is the desired number of PSWW modes, csi is the vector kx of the
% Floquet modes used to represent the PSWW modes, E, F are the Fourier
% transform of the square and 1/square refractive index profiles.
%
% kzTM is the vector of TM propagation constant, kTG2 is the vector of TM
% eigenvalues (mind the square!), VTM are the TM eigenvectors, ITE are
% computed from VTM. VTM, ITM have NModes rows and length(csi) columns
% (change of basis from Floquet mode basis to PSWW mode basis)
%
% Alberto Tibaldi, 04/03/2015
%------------------------------------------------------------------------

%-- Loading constants
f_LoadConstants

omega=k0.*Clight;

iF=inv(F);
iE=inv(E);

%-- Matrix of the eigenvalue problem (Lifeng Li formulation)
A=iF*(-diag(csi)*iE*diag(csi)+k0.^2.*eye(length(csi))); %-- Lifeng Li
% A=E*(-diag(csi)*iE*diag(csi)+k0.^2.*eye(length(csi))); %-- Moharam
% A=E*(-diag(csi)*F*diag(csi)+k0.^2.*eye(length(csi))); %-- "Classic"

[EigVc,EigVl]=eig(A);
EigVl=diag(EigVl);
%-- Sorting eigenvalues and eigenvectors
[BE, IN] = sort(real(EigVl));
IN = flipud(IN);
IN = IN(1:NModes);
EigVl = EigVl(IN);

indkT = find(angle(EigVl) > 0 & angle(EigVl) < 1e-5);
EigVl(indkT) = real(EigVl(indkT));
indkT = find(angle(EigVl) > (pi-1e-5) & angle(EigVl) < pi);
EigVl(indkT) = real(EigVl(indkT));
indkT = find(angle(EigVl) > 1e-5 & angle(EigVl) < (pi-1e-5));
if length(indkT) > 0 
	disp('WARNING!!! ');
	disp('kT2 CON PARTE IMMAGINARIA POSITIVA');
	disp('LA STRUTTURA E'' VISTA COME ATTIVA');
	disp(' ');
end

%-- Defining outputs
kTG2=EigVl;
kzTM=f_psqrt(kTG2-ky.^2);

ITM=EigVc(:,IN); %-- This is I'

matD=-j.*diag(csi);
NI=1./sqrt(diag(ITM.'*F*ITM)); %-- Orthonormalization is weighted!!!
ITM=ITM*diag(NI);

VTM=matD*ITM;
VTM=j/(k0/(Clight*mu0))*VTM;  %-- This is V', according to Tx line eqs.

return
