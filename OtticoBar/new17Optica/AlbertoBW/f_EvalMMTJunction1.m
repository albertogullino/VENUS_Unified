
function [S11,S21,S12,S22]=f_EvalMMTJunction1(LayerInfo,indLayer,Parameters,HalfSpaceInfo,ky)
%
% [S11,S21,S12,S22]=f_EvalMMTJunction1(lambda,LayerInfo,indLayer,Geometry,Parameters,HalfSpaceInfo,ky)
%
%------------------------------------------------------------------------
% GSM of the input junction
%------------------------------------------------------------------------
% This function computes, through a mode-matching technique (MMT), the
% generalized scattering matrix (GSM) of the junction between the left
% half-space and the first layer. "G" stands for "Guide", so for the PSWW,
% whereas "F" for Floquet, so for the half-space.
%
% Alberto Tibaldi, 04/03/2015
%------------------------------------------------------------------------

%-- Computation of the MMT projection matrices
[TeGF,ThFG]=f_EvalMMTProjMatJunction1(ky,Parameters,HalfSpaceInfo,LayerInfo,indLayer);

YF=diag([HalfSpaceInfo(1).YinfTE;HalfSpaceInfo(1).YinfTM]);
YG=diag(1./([LayerInfo(indLayer).ZG_TE;LayerInfo(indLayer).ZG_TM]));

ZF=diag(1./diag(YF));
ZG=diag(1./diag(YG));

%-- MMT
iW=inv(YF+ThFG*YG*TeGF);
I=eye(size(iW));
S11=sqrt(YF)*iW*(YF-ThFG*YG*TeGF)*sqrt(ZF);
S21=sqrt(YG)*TeGF*(I+iW*(YF-ThFG*YG*TeGF))*sqrt(ZF);
S12=2*sqrt(YF)*iW*ThFG*sqrt(YG);
S22=2*sqrt(YG)*TeGF*iW*ThFG*sqrt(YG)-eye(size(YG));

return