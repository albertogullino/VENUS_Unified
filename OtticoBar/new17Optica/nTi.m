function [n,k] = nTi( lambda )
%NCR( Lambda ) Returns the refractive index of Ti
%       Lambda is the wavelength in m.
% Curve fit to data from J.H. Weaver and H. P. R. Frederkse
% "Opticalproperties of metals and semiconductors"
% Valid in the interval 689 - 1240 nm.

if min(lambda) < 689e-9
        disp('Too short wavelength.');
        return
        elseif max(lambda) > 1.24e-6
        disp('Too large wavelength.');
%        return
        lambda=1.24e-6
end

e=6.6260755e-34*3e8./([1 1.10 1.20 1.30 1.40 1.50 1.60 1.70 1.80]*1.60217733e-19);
ntab=[3.62 3.47 3.35 3.28 3.17 2.98 2.74 2.54 2.36];
ktab=[3.52 3.40 3.30 3.25 3.28 3.32 3.30 3.23 3.11];

n = interp1(e,ntab,lambda,'spline');
k = interp1(e,ktab,lambda,'spline');

