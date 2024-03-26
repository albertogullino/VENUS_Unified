function nindex=nAlGaAs(lambda,x)
%NALGAAS(LAMBDA,X) Returns the refractive index of AlGaAs
%       LAMBDA is the wavelength in m (can be a vector) and
%       X is the compositional Al content.
%       NALGAAS(LAMBDA,X) returns an array with the same size as LAMBDA
%
%       Original model from M.A. Afromovitz, Solid State Comm, 15, 1974
%       Model extended with the range 0.8µm - egamma from
%       S.W. Corzine, Ph.D. Dissertation

% iAlGaAs_model = 1;  % 1 -> our model, 2 -> TRUMPF model
iAlGaAs_model = 2;  % 1 -> our model, 2 -> TRUMPF model

if iAlGaAs_model == 1
    nindex=nAlGaAs_Corzine(lambda,x);
else
    nindex=nAlGaAs_gehrsitz_WL(lambda,x);
end

end