function [n,k] = nPt( lambda )
%NCR( Lambda ) Returns the refractive index of Pt
%       Lambda is the wavelength in m.
% Curve fit to data from J.H. Weaver and H. P. R. Frederkse
% "Opticalproperties of metals and semiconductors"
% Valid in the interval 689 - 1240 nm.

if min(lambda) < 689e-9
        disp('Too short wavelength.');
        return
        elseif max(lambda) > 1.24e-6
        disp('Too large wavelength.');
        lambda=1.24e-6;
        %return
end

e=6.6260755e-34*3e8./([1 1.10 1.20 1.30 1.40 1.50 1.60 1.70 1.80]*1.60217733e-19);
ntab=[4.25 3.86 3.55 3.29 3.10 2.92 2.76 2.63 2.51];
ktab=[6.62 6.24 5.92 5.61 5.32 5.07 4.84 4.64 4.43];

n = interp1(e,ntab,lambda,'spline');
k = interp1(e,ktab,lambda,'spline');

