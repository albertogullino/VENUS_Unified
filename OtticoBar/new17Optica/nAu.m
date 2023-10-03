function [n,k] = nAu( lambda )
%NAU( Lambda ) Returns the refractive index of Au
%       Lambda is the wavelength in m.
% Curve fit to data from Palik "Handbook of Optical
% Constants of Solids

n = 8.39334e-2 - lambda*8.09271e4 + lambda.^2*2.52477e11;
k = -3.05320 + lambda*1.14275e7 - lambda.^2*1.52222e12;

%k =( -3.05320 + lambda*1.14275e7 - lambda.^2*1.52222e12) /10;
%n=0.1;
%k=10;
%' attenzione nAu modificato ', keyboard
