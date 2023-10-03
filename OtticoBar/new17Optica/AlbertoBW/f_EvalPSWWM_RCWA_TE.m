function [kzTE,kTG2,VTE,ITE]=f_EvalPSWWM_RCWA_TE(k0,NModes,csi,ky,E)
%
% [kzTE,kTG2,VTE,ITE]=f_EvalPSWWM_RCWA_TE(k0,NModes,csi,ky,E)
%
%------------------------------------------------------------------------
% Synthesis of PSWW TE modes
%------------------------------------------------------------------------
% This function calculates the propagation constants (along z) and the
% change of basis matrices from the generic Floquet modes basis to the PSWW
% modes, as the solution of an eigenvalue problem. This eigenvalue problem
% is regular, therefore no special formulation is implemented.
%
% NModes is the desired number of PSWW modes, csi is the vector kx of the
% Floquet modes used to represent the PSWW modes, E is the Fourier
% transform of the square refractive index profile.
%
% kzTE is the vector of TE propagation constant, kTG2 is the vector of TE
% eigenvalues (mind the square!), VTE are the TE eigenvectors, ITE are
% computed from VTE. VTE, ITE have NModes rows and length(csi) columns
% (change of basis from Floquet mode basis to PSWW mode basis)
%
% Alberto Tibaldi, 04/03/2015
%------------------------------------------------------------------------

%-- Loading constants
f_LoadConstants

omega=k0.*Clight;

%-- Matrix of the eigenvalue problem
A=k0.^2.*E-diag(csi.^2);

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
kzTE=f_psqrt(kTG2-ky.^2);

VTE=EigVc(:,IN); %-- This is V''

matD=-j.*diag(csi);
ITE=matD*VTE;
ITE=j/(omega.*mu0)*ITE; %-- This is I'', according to Tx line equations

return
