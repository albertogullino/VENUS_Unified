function [n,k] = nCr( lambda )
%NCR( Lambda ) Returns the refractive index of Cr
%	Lambda is the wavelength in m.
% Curve fit to data from J.H. Weaver and H. P. R. Frederkse
% "Opticalproperties of metals and semiconductors"
% Valid in the interval 700 - 1240 nm.

if min(lambda) < 700e-9
	disp('Too short wavelength.');
	return
	elseif max(lambda) > 1.24e-6
	disp('Too large wavelength.');
	return
end

e=6.6260755e-34*3e8./([1 1.12 1.24 1.36 1.46 1.77]*1.60217733e-19);
ntab=[4.47 4.53 4.50 4.42 4.31 3.84];
ktab=[4.43 4.31 4.28 4.30 4.32 4.37];

n = interp1(e,ntab,lambda,'spline');
k = interp1(e,ktab,lambda,'spline');

% np = [4.496792644642644e+18    -1.760878059828245e+13
2.231578333243547e+07...
%     -4.698772806816330e+00];
% kp = [3.263225001126466e+18    -7.940349666179611e+12
6.064592611072156e+06...
%      2.896867746659526e+00];
% n = np(4) + lambda*np(3) + lambda.^2*np(2) + lambda.^3*np(1);
% k = kp(4) + lambda*kp(3) + lambda.^2*kp(2) + lambda.^3*kp(1);
