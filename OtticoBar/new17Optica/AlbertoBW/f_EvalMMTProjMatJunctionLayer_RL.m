
function [TeRL,ThLR]=f_EvalMMTProjMatJunctionLayer_RL(lambda,LayerInfo,indLayer,Parameters,ky)
%
% [TeRL,ThLR]=f_EvalMMTProjMatJunctionLayer(lambda,LayerInfo,indLayer,Parameters,Geometry,theta,phi)
%
%------------------------------------------------------------------------
% MMT Projection matrices junction among two layers
%------------------------------------------------------------------------
% This function returns the projection matrices for the MMT aimed at
% characterizing the junction between two PSWWs.
% This routine is related to the (R,L) formulation.
%
% Alberto Tibaldi, 04/03/2015
%------------------------------------------------------------------------

f_LoadConstants
k0=2*pi/lambda;
omega=k0*Clight;

d=LayerInfo(indLayer).d1+LayerInfo(indLayer).d2;
d1=LayerInfo(indLayer-1).d1;

ITERH=LayerInfo(indLayer).ITEH(1:LayerInfo(indLayer).NModes,1:Parameters.NHarmonicsTE);
VTERH=LayerInfo(indLayer).VTEH(1:LayerInfo(indLayer).NModes,1:Parameters.NHarmonicsTE);
ITMRH=LayerInfo(indLayer).ITMH(1:LayerInfo(indLayer).NModes,1:Parameters.NHarmonicsTM);
ITER=LayerInfo(indLayer).ITE(1:LayerInfo(indLayer).NModes,1:Parameters.NHarmonicsTE);
VTER=LayerInfo(indLayer).VTE(1:LayerInfo(indLayer).NModes,1:Parameters.NHarmonicsTE);
ITMR=LayerInfo(indLayer).ITM(1:LayerInfo(indLayer).NModes,1:Parameters.NHarmonicsTM);

VTELH=LayerInfo(indLayer-1).VTEH(1:LayerInfo(indLayer-1).NModes,1:Parameters.NHarmonicsTE);
ITMLH=LayerInfo(indLayer-1).ITMH(1:LayerInfo(indLayer-1).NModes,1:Parameters.NHarmonicsTM);
VTMLH=LayerInfo(indLayer-1).VTMH(1:LayerInfo(indLayer-1).NModes,1:Parameters.NHarmonicsTM);
VTEL=LayerInfo(indLayer-1).VTE(1:LayerInfo(indLayer-1).NModes,1:Parameters.NHarmonicsTE);
ITML=LayerInfo(indLayer-1).ITM(1:LayerInfo(indLayer-1).NModes,1:Parameters.NHarmonicsTM);
VTML=LayerInfo(indLayer-1).VTM(1:LayerInfo(indLayer-1).NModes,1:Parameters.NHarmonicsTM);

csiLTE=LayerInfo(indLayer-1).csiHarmTE;
csiRTE=LayerInfo(indLayer).csiHarmTE;
csiLTM=LayerInfo(indLayer-1).csiHarmTM;
csiRTM=LayerInfo(indLayer).csiHarmTM;
n1L=LayerInfo(indLayer-1).n1;
n2L=LayerInfo(indLayer-1).n2;

%-- Raised cosine profile parameters (temp)

[NHarmMax,ind]=max([Parameters.NHarmonicsTE,Parameters.NHarmonicsTM]);
if ind==1
    FLTETM=LayerInfo(indLayer-1).FTE(:,1:Parameters.NHarmonicsTM);
    FLTMTE=LayerInfo(indLayer-1).FTE(1:Parameters.NHarmonicsTM,:);
    FLTMTM=LayerInfo(indLayer-1).FTE(1:Parameters.NHarmonicsTM,1:Parameters.NHarmonicsTM);
elseif ind==2
    FLTETM=LayerInfo(indLayer-1).FTM(1:Parameters.NHarmonicsTE,:);
    FLTMTE=LayerInfo(indLayer-1).FTM(:,1:Parameters.NHarmonicsTE);
    FLTMTM=LayerInfo(indLayer-1).FTM;
end

TeRTELTE=VTERH*(VTEL).';
TeRTELTMA=-ky.*omega.*eps0.*VTERH*(diag(1./LayerInfo(indLayer-1).kTTM.^2)*VTML*FLTETM.').';
TeRTELTMB=-omega.*mu0.*ky.*(diag(1./LayerInfo(indLayer).kTTE.^2)*ITERH)*(ITML*FLTETM.').';
TeRTELTM=TeRTELTMA+TeRTELTMB;
TeRTMLTE=zeros(LayerInfo(indLayer).NModes,LayerInfo(indLayer-1).NModes);
TeRTMLTM=ITMRH*(ITML*FLTMTM.').';

TeRL=[TeRTELTE,TeRTELTM;
    TeRTMLTE,TeRTMLTM];

ThLTMRTE1=omega.*mu0.*ky.*((ITMLH*FLTMTE)*(diag(1./LayerInfo(indLayer).kTTE.^2)*ITER).');
ThLTMRTE2=omega.*eps0.*ky.*((diag(1./LayerInfo(indLayer-1).kTTM.^2)*VTMLH*FLTMTE)*(VTER).');
ThLTMRTE=ThLTMRTE1+ThLTMRTE2;
ThLTERTE=VTELH*VTER.';
ThLTERTM=zeros(LayerInfo(indLayer-1).NModes,LayerInfo(indLayer).NModes);
ThLTMRTM=(ITMLH*FLTMTM)*ITMR.';

ThLR=[ThLTERTE,ThLTERTM;
    ThLTMRTE,ThLTMRTM];

return