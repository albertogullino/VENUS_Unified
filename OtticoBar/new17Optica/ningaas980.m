function nOut = nInGaAs980(lambda, x)
%NINGAAS980(Lambda,x) Returns the refractive index of In(x)Ga(1-x)As
%	Lambda is the wavelength in m, should be between 0.8856µm and
%1.033µm, and
%	x is the compositional In content, can be a vector.

n_InAs_08856 = 3.696;
n_InAs_1033 = 3.613;
n_InAs = n_InAs_08856 + (n_InAs_1033 - n_InAs_08856)/(1.033 - .8856)...
			*(lambda*1e6 - .8856);
n_GaAs = nAlGaAs( lambda, 0 );

nOut = x*n_InAs + (1-x)*n_GaAs;
