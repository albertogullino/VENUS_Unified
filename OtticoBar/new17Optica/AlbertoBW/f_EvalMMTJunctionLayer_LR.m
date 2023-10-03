
function [S11,S21,S12,S22]=f_EvalMMTJunctionLayer_LR(lambda,LayerInfo,indLayer,Parameters,ky)
%
% [TeGF,ThFG]=f_EvalMMTProjMatJunction1(ky,Parameters,HalfSpaceInfo,LayerInfo,indLayer)
%
%------------------------------------------------------------------------
% MMT Projection matrices junction among two layers
%------------------------------------------------------------------------
% This function returns the projection matrix for the MMT aimed at
% characterizing the junction between two inner layers. "L" stands for
% "Left" tooth, wheras "R" stands for "Right".
% This routine is related to the (L,R) formulation.
%
% Alberto Tibaldi, 04/03/2015
%------------------------------------------------------------------------

%-- Computation of the MMT projection matrices
[TeLR,ThRL]=f_EvalMMTProjMatJunctionLayer_LR(lambda,LayerInfo,indLayer,Parameters,ky);

YL=diag(1./([LayerInfo(indLayer-1).ZG_TE;LayerInfo(indLayer-1).ZG_TM]));
YR=diag(1./([LayerInfo(indLayer).ZG_TE;LayerInfo(indLayer).ZG_TM]));

ZL=diag(1./diag(YL));
ZR=diag(1./diag(YR));

%-- MMT
iW=inv(YR+ThRL*YL*TeLR);
IL=eye(size(YL));
IR=eye(size(YR));
S11=2*sqrt(YL)*TeLR*iW*ThRL*sqrt(YL)-IL;
S21=2*sqrt(YR)*iW*ThRL*sqrt(YL);
S12=sqrt(YL)*TeLR*(IR+iW*(YR-ThRL*YL*TeLR))*sqrt(ZR);
S22=sqrt(YR)*iW*(YR-ThRL*YL*TeLR)*sqrt(ZR);

return