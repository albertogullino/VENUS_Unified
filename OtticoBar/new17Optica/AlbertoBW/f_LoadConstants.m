
%------------------------------------------------------------------------
% This file contains the universal constants, derived starting from the
% definition of Clight and of mu0. By calling it, these quantities are
% stored in proper variables.
%
% Alberto Tibaldi, 18/02/2015
%------------------------------------------------------------------------

Clight=299.792458; % speed of light (um/ps)
mu0=4*pi/10; % vacuum magnetic permeability (H/um)
Z0=mu0*Clight; % free-space impedance (ohm)
eps0=1./(mu0.*Clight.^2);