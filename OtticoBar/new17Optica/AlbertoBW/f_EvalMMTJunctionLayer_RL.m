
function [S11,S21,S12,S22]=f_EvalMMTJunctionLayer_RL(lambda,LayerInfo,indLayer,Parameters,ky)
%
% [TeGF,ThFG]=f_EvalMMTProjMatJunction1(ky,Parameters,HalfSpaceInfo,LayerInfo,indLayer)
%
%------------------------------------------------------------------------
% MMT Projection matrices junction among two layers
%------------------------------------------------------------------------
% This function returns the projection matrix for the MMT aimed at
% characterizing the junction between two inner layers. "L" stands for
% "Left" tooth, wheras "R" stands for "Right".
% This routine is related to the (R,L) formulation.
%
% Alberto Tibaldi, 04/03/2015
%------------------------------------------------------------------------

%-- Computation of the MMT projection matrices
[TeRL,ThLR]=f_EvalMMTProjMatJunctionLayer_RL(lambda,LayerInfo,indLayer,Parameters,ky);

YL=diag(1./([LayerInfo(indLayer-1).ZG_TE;LayerInfo(indLayer-1).ZG_TM]));
YR=diag(1./([LayerInfo(indLayer).ZG_TE;LayerInfo(indLayer).ZG_TM]));

ZL=diag(1./diag(YL));
ZR=diag(1./diag(YR));

%-- MMT
iW=inv(YL+ThLR*YR*TeRL);
I=eye(size(iW));
S11=sqrt(YL)*iW*(YL-ThLR*YR*TeRL)*sqrt(ZL);
S21=sqrt(YR)*TeRL*(I+iW*(YL-ThLR*YR*TeRL))*sqrt(ZL);
S12=2*sqrt(YL)*iW*ThLR*sqrt(YR);
S22=2*sqrt(YR)*TeRL*iW*ThLR*sqrt(YR)-eye(size(YR));

return