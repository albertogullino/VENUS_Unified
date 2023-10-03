
function [TeLR,ThRL]=f_EvalMMTProjMatJunctionLayer_LR(lambda,LayerInfo,indLayer,Parameters,ky)
%
% [TeRL,ThLR]=f_EvalMMTProjMatJunctionLayer(lambda,LayerInfo,indLayer,Parameters,Geometry,theta,phi)
%
%------------------------------------------------------------------------
% MMT Projection matrices junction among two layers
%------------------------------------------------------------------------
% This function returns the projection matrices for the MMT aimed at
% characterizing the junction between two PSWWs.
%
% Alberto Tibaldi, 04/03/2015
%------------------------------------------------------------------------

f_LoadConstants
k0=2*pi/lambda;
omega=k0*Clight;

d=LayerInfo(indLayer).d1+LayerInfo(indLayer).d2;
d1=LayerInfo(indLayer).d1;

ITELH=LayerInfo(indLayer-1).ITEH(1:LayerInfo(indLayer-1).NModes,1:Parameters.NHarmonicsTE);
VTERH=LayerInfo(indLayer).VTEH(1:LayerInfo(indLayer).NModes,1:Parameters.NHarmonicsTE);
ITMRH=LayerInfo(indLayer).ITMH(1:LayerInfo(indLayer).NModes,1:Parameters.NHarmonicsTM);
ITEL=LayerInfo(indLayer-1).ITE(1:LayerInfo(indLayer-1).NModes,1:Parameters.NHarmonicsTE);
VTER=LayerInfo(indLayer).VTE(1:LayerInfo(indLayer).NModes,1:Parameters.NHarmonicsTE);
ITMR=LayerInfo(indLayer).ITM(1:LayerInfo(indLayer).NModes,1:Parameters.NHarmonicsTM);

VTELH=LayerInfo(indLayer-1).VTEH(1:LayerInfo(indLayer-1).NModes,1:Parameters.NHarmonicsTE);
ITMLH=LayerInfo(indLayer-1).ITMH(1:LayerInfo(indLayer-1).NModes,1:Parameters.NHarmonicsTM);
VTMRH=LayerInfo(indLayer).VTMH(1:LayerInfo(indLayer).NModes,1:Parameters.NHarmonicsTM);
VTEL=LayerInfo(indLayer-1).VTE(1:LayerInfo(indLayer-1).NModes,1:Parameters.NHarmonicsTE);
ITML=LayerInfo(indLayer-1).ITM(1:LayerInfo(indLayer-1).NModes,1:Parameters.NHarmonicsTM);
VTMR=LayerInfo(indLayer).VTM(1:LayerInfo(indLayer).NModes,1:Parameters.NHarmonicsTM);

[NHarmMax,ind]=max([Parameters.NHarmonicsTE,Parameters.NHarmonicsTM]);
if ind==1
    FRTETM=LayerInfo(indLayer).FTE(:,1:Parameters.NHarmonicsTM);
    FRTMTE=LayerInfo(indLayer).FTE(1:Parameters.NHarmonicsTM,:);
    FRTMTM=LayerInfo(indLayer).FTE(1:Parameters.NHarmonicsTM,1:Parameters.NHarmonicsTM);
elseif ind==2
    FRTETM=LayerInfo(indLayer).FTM(1:Parameters.NHarmonicsTE,:);
    FRTMTE=LayerInfo(indLayer).FTM(:,1:Parameters.NHarmonicsTE);
    FRTMTM=LayerInfo(indLayer).FTM;
end

TeLTERTE=VTELH*(VTER).';
TeLTERTMA=-ky.*omega.*eps0.*VTELH*(diag(1./LayerInfo(indLayer).kTTM.^2)*VTMR*FRTETM.').';
TeLTERTMB=-omega.*mu0.*ky.*(diag(1./LayerInfo(indLayer-1).kTTE.^2)*ITELH)*(ITMR*FRTETM.').';
TeLTERTM=TeLTERTMA+TeLTERTMB;
TeLTMRTE=zeros(LayerInfo(indLayer-1).NModes,LayerInfo(indLayer).NModes);
TeLTMRTM=ITMLH*(ITMR*FRTMTM.').';

TeLR=[TeLTERTE,TeLTERTM;
    TeLTMRTE,TeLTMRTM];

ThRTMLTE1=omega.*mu0.*ky.*((ITMRH*FRTMTE)*(diag(1./LayerInfo(indLayer-1).kTTE.^2)*ITEL).');
ThRTMLTE2=omega.*eps0.*ky.*((diag(1./LayerInfo(indLayer).kTTM.^2)*VTMRH*FRTMTE)*(VTEL).');
ThRTMLTE=ThRTMLTE1+ThRTMLTE2;
ThRTELTE=VTERH*VTEL.';
ThRTELTM=zeros(LayerInfo(indLayer).NModes,LayerInfo(indLayer-1).NModes);
ThRTMLTM=(ITMRH*FRTMTM)*ITML.';

ThRL=[ThRTELTE,ThRTELTM;
    ThRTMLTE,ThRTMLTM];

return