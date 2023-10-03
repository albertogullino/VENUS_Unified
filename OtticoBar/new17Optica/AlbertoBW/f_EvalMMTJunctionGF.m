
function [S11,S21,S12,S22]=f_EvalMMTJunctionGF(LayerInfo,indLayer,Parameters,HalfSpaceInfo,ky)
%
% [S11,S21,S12,S22]=f_EvalMMTJunctionGF(LayerInfo,indLayer,Parameters,HalfSpaceInfo,ky)
%
%------------------------------------------------------------------------
% GSM from grating (G) to homogeneous waveguide (F)
%------------------------------------------------------------------------
%
% Alberto Tibaldi, 04/03/2015
%------------------------------------------------------------------------

[TeLR,ThRL]=f_EvalMMTProjMatJunctionEnd(ky,Parameters,HalfSpaceInfo,LayerInfo,indLayer-1);
%-- In questo caso YF è quella di destra, YG quella di sinistra, stai
%attento!!!!

YR=diag(1./([LayerInfo(indLayer).ZG_TE;LayerInfo(indLayer).ZG_TM]));
YL=diag(1./([LayerInfo(indLayer-1).ZG_TE;LayerInfo(indLayer-1).ZG_TM]));

ZR=diag(1./diag(YR));
ZL=diag(1./diag(YL));

iW=inv(YR+ThRL*YL*TeLR);
I=eye(size(iW));
S22=sqrt(YR)*iW*(YR-ThRL*YL*TeLR)*sqrt(ZR);
S12=sqrt(YL)*TeLR*(I+iW*(YR-ThRL*YL*TeLR))*sqrt(ZR);
S21=2*sqrt(YR)*iW*ThRL*sqrt(YL);
S11=2*sqrt(YL)*TeLR*iW*ThRL*sqrt(YL)-eye(size(YL));