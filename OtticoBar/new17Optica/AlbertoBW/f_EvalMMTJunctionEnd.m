
function [S11,S21,S12,S22]=f_EvalMMTJunctionEnd(LayerInfo,indLayer,Parameters,HalfSpaceInfo,ky)
%
% [S11,S21,S12,S22]=f_EvalMMTJunctionEnd(LayerInfo,indLayer,Parameters,HalfSpaceInfo,ky)
%
%------------------------------------------------------------------------
% GSM of the ending junction
%------------------------------------------------------------------------
% This function computes, through a mode-matching technique (MMT), the
% generalized scattering matrix (GSM) of the junction between the right
% half-space and the last layer. "G" stands for "Guide", so for the PSWW,
% whereas "F" for Floquet, so for the half-space.
%
% Alberto Tibaldi, 04/03/2015
%------------------------------------------------------------------------

[TeGF,ThFG]=f_EvalMMTProjMatJunctionEnd(ky,Parameters,HalfSpaceInfo,LayerInfo,indLayer);

YF=diag([HalfSpaceInfo(2).YinfTE;HalfSpaceInfo(2).YinfTM]);
YG=diag(1./([LayerInfo(indLayer).ZG_TE;LayerInfo(indLayer).ZG_TM]));

ZF=diag(1./diag(YF));
ZG=diag(1./diag(YG));

iW=inv(YF+ThFG*YG*TeGF);
I=eye(size(iW));
S22=sqrt(YF)*iW*(YF-ThFG*YG*TeGF)*sqrt(ZF);
S12=sqrt(YG)*TeGF*(I+iW*(YF-ThFG*YG*TeGF))*sqrt(ZF);
S21=2*sqrt(YF)*iW*ThFG*sqrt(YG);
S11=2*sqrt(YG)*TeGF*iW*ThFG*sqrt(YG)-eye(size(YG));