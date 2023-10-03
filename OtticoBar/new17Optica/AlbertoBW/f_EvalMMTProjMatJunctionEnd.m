
function [TeGF,ThFG]=f_EvalMMTProjMatJunctionEnd(ky,Parameters,HalfSpaceInfo,LayerInfo,indLayer)
%
% [TeGF,ThFG]=f_EvalMMTProjMatJunctionEnd(ky,Parameters,HalfSpaceInfo,LayerInfo,indLayer)
%
%------------------------------------------------------------------------
% MMT Projection matrices last junction
%------------------------------------------------------------------------
% This function returns the projection matrices for the MMT aimed at
% characterizing the junction between the right half-space and the last
% layer. The integrals are computed analytically.
%
% Alberto Tibaldi, 04/03/2015
%------------------------------------------------------------------------

ITEH=LayerInfo(indLayer).ITEH(1:LayerInfo(indLayer).NModes,1:Parameters.NModes);
VTEH=LayerInfo(indLayer).VTEH(1:LayerInfo(indLayer).NModes,1:Parameters.NModes);
ITMH=LayerInfo(indLayer).ITMH(1:LayerInfo(indLayer).NModes,1:Parameters.NModes);
ITE=LayerInfo(indLayer).ITE(1:LayerInfo(indLayer).NModes,1:Parameters.NModes);
VTE=LayerInfo(indLayer).VTE(1:LayerInfo(indLayer).NModes,1:Parameters.NModes);
ITM=LayerInfo(indLayer).ITM(1:LayerInfo(indLayer).NModes,1:Parameters.NModes);

kx_on_kTF=HalfSpaceInfo(2).kx./HalfSpaceInfo(2).kTF.';
ky_on_kTF=ky./HalfSpaceInfo(2).kTF;
ind0=find(HalfSpaceInfo(2).kTF==0);
kx_on_kTF(ind0)=1;
ky_on_kTF(ind0)=0;

% TGTEFTEe=-j.*ky.*diag(LayerInfo(indLayer).ZG_TE./LayerInfo(indLayer).kzTE)*ITEH*diag(ky_on_kTF)-j.*VTEH*diag(kx_on_kTF);
TGTEFTEe=-ky.*diag(LayerInfo(indLayer).ZG_TE./LayerInfo(indLayer).kzTE)*ITEH*diag(ky_on_kTF)-VTEH*diag(kx_on_kTF);
TGTEFTMe=-j.*ky.*diag(LayerInfo(indLayer).ZG_TE./LayerInfo(indLayer).kzTE)*ITEH*diag(kx_on_kTF)+j.*VTEH*diag(ky_on_kTF);
TGTMFTEe=j.*ITMH*diag(ky_on_kTF);
% TGTMFTMe=j.*ITMH*diag(kx_on_kTF);
TGTMFTMe=ITMH*diag(kx_on_kTF);

% TFTEGTEh=j.*(VTE*diag(kx_on_kTF)).'-j.*ky.*(diag(LayerInfo(indLayer).ZG_TE./LayerInfo(indLayer).kzTE)*ITE*diag(ky_on_kTF)).';
TFTEGTEh=-(VTE*diag(kx_on_kTF)).'+ky.*(diag(LayerInfo(indLayer).ZG_TE./LayerInfo(indLayer).kzTE)*ITE*diag(ky_on_kTF)).';
TFTMGTEh=-j.*(VTE*diag(ky_on_kTF)).'-j.*ky.*(diag(LayerInfo(indLayer).ZG_TE./LayerInfo(indLayer).kzTE)*ITE*diag(kx_on_kTF)).';
TFTEGTMh=-j.*(ITM*diag(ky_on_kTF)).';
% TFTMGTMh=-j.*(ITM*diag(kx_on_kTF)).';
TFTMGTMh=(ITM*diag(kx_on_kTF)).';

ThFG=[(TFTEGTEh),(TFTEGTMh);
    (TFTMGTEh),(TFTMGTMh)];

TeGF=[(TGTEFTEe),(TGTEFTMe);
    (TGTMFTEe),(TGTMFTMe)];

return