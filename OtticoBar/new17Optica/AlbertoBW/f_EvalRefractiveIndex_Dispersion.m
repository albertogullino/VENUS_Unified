function [n]=f_EvalRefractiveIndex(lambda,material)
%
% function [PatchInfo] = f_GetPCMatrices(PatchInfo,PCInfo,Parameters,Geometry,nodes,weights,indf)
%
%------------------------------------------------------------------------
% Evaluation of the refractive index
%------------------------------------------------------------------------
% This function returns the refractive index of a certain material, for a
% chosen wavelength lambda.
%
% lambda is the free-space wavelength
%
% material is a string containing the material.
%
% Available materials are: 'Si', 'SiO2', 'Pyrex', 'SiNonDisp',
% 'SiO2NonDisp', 'PyrexNonDisp', 'Vacuum'
%
% Alberto Tibaldi, 18/02/2015
%------------------------------------------------------------------------

if strcmp(material,'Si')==1
    %-- Amorphous silicon
    %-- Approximation between Si and amorphous Si, + experiments
    n = 1.05.*sqrt(1 + 10.6684293*lambda.^2./(lambda.^2-0.301516485^2) + 0.003043475*lambda.^2./(lambda.^2-1.13475115^2) + 1.54133408*lambda.^2./(lambda.^2-1104.0^2));
elseif strcmp(material,'SiO2')==1
    %-- Silicon dioxide
    %-- http://refractiveindex.info/?group=CRYSTALS&material=SiO2
    n = sqrt( 1 + 0.663044.*lambda.^2./(lambda.^2-0.060^2) + 0.517852.*lambda.^2./(lambda.^2-0.106^2) + 0.175912.*lambda.^2./(lambda.^2-0.119^2) + 0.565380.*lambda.^2./(lambda.^2-8.844^2) + 1.675299.*lambda.^2./(lambda.^2-20.742^2) );
elseif strcmp(material,'Pyrex')==1
    %-- Borosilicate
    %-- http://refractiveindex.info/?shelf=glass&book=SCHOTT-multipurpose&page=BOROFLOAT33.yml
    %-- extrapolation in infrared, with spline
    n = spline([0.4358, 0.4799, 0.5461, 0.5871, 0.5893, 0.6438, 0.6563, 1, 2, 3, 4],[1.48015, 1.47676, 1.47311, 1.47140, 1.47133, 1.46953, 1.46916, 1.46, 1.445, 1.437, 1.433], lambda);
elseif strcmp(material,'SiNonDisp')==1
    n=3.43;
elseif strcmp(material,'SiO2NonDisp')==1
    n=1.41;
elseif strcmp(material,'PyrexNonDisp')==1
    n=1.47;
elseif strcmp(material,'Vacuum')==1
    n=1;
else
    disp('Material not available; returning n=0')
    n=0;
    return
end

return