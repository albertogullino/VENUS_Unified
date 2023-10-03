function [LayerInfo]=f_EvalPSWW_Modes(lambda,indLayer,Parameters,Geometry,LayerInfo,Options,kx0,ky)
%
% [LayerInfo]=f_EvalPSWW_Modes(lambda,indLayer,Parameters,Geometry,LayerInfo,Options,KBd,ky)
%
%------------------------------------------------------------------------
% Synthesis of PSWW modes
%------------------------------------------------------------------------
% This function is used to compute the phase-shift wall waveguide modes
% (PSWW modes) for the LayerInfo(indLayer) structure. The refractive index
% profile in the waveguide is a raised-cosine. The PSWW modes are
% calculating according to the rigorous coupled wave analysis (RCWA), that
% consists of representing the modes of the structure in terms of Floquet
% modes.
%
% Alberto Tibaldi, 04/03/2015
%------------------------------------------------------------------------

%-- Loading constants
f_LoadConstants

k0=2*pi/lambda;
omega=k0*Clight;
eps0=k0/(Clight*omega*mu0);

d=LayerInfo(indLayer).d1+LayerInfo(indLayer).d2;
d1=LayerInfo(indLayer).d1;

%-- Raised cosine profile parameters (temp)
alpha1=LayerInfo(indLayer).RollOff;
a1=d1/2;
x1=LayerInfo(indLayer).Displacement;

%-- Defining the Floquet modes basis used to represent TE PSWW modes
vnTE_temp = 1:(Parameters.NHarmonicsTE-1)/2;
vnTE(2*vnTE_temp) = -vnTE_temp;
vnTE(2*vnTE_temp+1) = vnTE_temp;

csiHarmTE=kx0 + 2*pi*vnTE./d; %-- kx of the expansion Floquet modes
LayerInfo(indLayer).vnTE=vnTE;
LayerInfo(indLayer).csiHarmTE=csiHarmTE;

%-- Defining the Floquet modes basis used to represent TM PSWW modes
vnTM_temp = 1:(Parameters.NHarmonicsTM-1)/2;
vnTM(2*vnTM_temp) = -vnTM_temp;
vnTM(2*vnTM_temp+1) = vnTM_temp;

csiHarmTM=kx0 + 2*pi*vnTM./d;
LayerInfo(indLayer).vnTM=vnTM;
LayerInfo(indLayer).csiHarmTM=csiHarmTM;

%-- Calculating Fourier/Floquet expansion of n^2(x) (E), 1/n^2(x) (F)

[NHarmMax,ind]=max([Parameters.NHarmonicsTE,Parameters.NHarmonicsTM]);
if ind==1
    if LayerInfo(indLayer).RollOff==0
        [ETE,FTE]=f_EvalRefractiveIndexExpansions(Parameters.NHarmonicsTE,csiHarmTE,d,alpha1,a1,x1,LayerInfo(indLayer).n1,LayerInfo(indLayer).n2);
    else
        N=Parameters.NFFT/2;
        vN=(-N:N-1); %-- Best choice, to maintain the matlab fft() order!xvet=vN*d/NFFT; %-- Sampling points in natural (spatial) domain
        xvet=vN/Parameters.NFFT*d;
        kxvet=vN/d; %-- Sampling points in the FFT reciprocal lattice (k)
        ExactProfile=LayerInfo(indLayer).n2^2+(LayerInfo(indLayer).n1^2-LayerInfo(indLayer).n2^2)*f_cosrialz(xvet,x1,d1,alpha1);
        [ETE]=f_EvalProfileFFT(ExactProfile,vN,Parameters.NHarmonicsTE,d);
        [FTE]=f_EvalProfileFFT(1./ExactProfile,vN,Parameters.NHarmonicsTE,d);
    end
    ETM=ETE(1:Parameters.NHarmonicsTM,1:Parameters.NHarmonicsTM);
    FTM=FTE(1:Parameters.NHarmonicsTM,1:Parameters.NHarmonicsTM);
    
elseif ind==2
    if LayerInfo(indLayer).RollOff==0
        [ETM,FTM]=f_EvalRefractiveIndexExpansions(Parameters.NHarmonicsTM,csiHarmTM,d,alpha1,a1,x1,LayerInfo(indLayer).n1,LayerInfo(indLayer).n2);
    else
        N=Parameters.NFFT/2;
        vN=(-N:N-1); %-- Best choice, to maintain the matlab fft() order!xvet=vN*d/NFFT; %-- Sampling points in natural (spatial) domain
        xvet=vN/Parameters.NFFT*d;
        ExactProfile=LayerInfo(indLayer).n2^2+(LayerInfo(indLayer).n1^2-LayerInfo(indLayer).n2^2)*f_cosrialz(xvet,x1,d1,alpha1);
        [ETM]=f_EvalProfileFFT(ExactProfile,vN,Parameters.NHarmonicsTM,d);
        [FTM]=f_EvalProfileFFT(1./ExactProfile,vN,Parameters.NHarmonicsTM,d);
    end
    ETE=ETM(1:Parameters.NHarmonicsTE,1:Parameters.NHarmonicsTE);
    FTE=FTM(1:Parameters.NHarmonicsTE,1:Parameters.NHarmonicsTE);    
end

LayerInfo(indLayer).ETE=ETE;
LayerInfo(indLayer).ETM=ETM;
LayerInfo(indLayer).FTE=FTE;
LayerInfo(indLayer).FTM=FTM;

if Options.RebuildProfile==1 %-- Drawing re-built profile

    x=linspace(-d/2,d/2,1001);
    %-- To draw the re-built profile, use Fourier, NOT Floquet, expansion
    FSMat=exp(-j*x.'*(csiHarmTE-csiHarmTE(1))); 
    figure(11),clf
    axis on
    grid on
    hold on
    RebuiltProfile=FSMat*ETE(:,1);
    plot(x,real(RebuiltProfile),'k')
    ExactProfile=LayerInfo(indLayer).n2^2+(LayerInfo(indLayer).n1^2-LayerInfo(indLayer).n2^2)*f_cosrialz(x,x1,d1,alpha1); 
    plot(x,ExactProfile)
    
end

%-- Determination of PSWW TE modes
[kzTE,kTG_TE2,VTE,ITE]=f_EvalPSWWM_RCWA_TE(k0,LayerInfo(indLayer).NModes,csiHarmTE,ky,ETE);

if imag(Geometry.n_in)==0 && imag(LayerInfo(indLayer).n1)==0 && imag(LayerInfo(indLayer).n2)==0
    %-- Lossless case
    VTEH=conj(VTE);
    ITEH=-conj(ITE); %-- "-" sign because of (j*xi*...) instead of (-j*xi*...)
else
    if LayerInfo(indLayer).RollOff==0
        [ETEH,FTEH]=f_EvalRefractiveIndexExpansions(Parameters.NHarmonicsTE,-csiHarmTE,d,alpha1,a1,x1,LayerInfo(indLayer).n1,LayerInfo(indLayer).n2);
    else
        N=Parameters.NFFT/2;
        vN=(-N:N-1); %-- Best choice, to maintain the matlab fft() order!xvet=vN*d/NFFT; %-- Sampling points in natural (spatial) domain
        xvet=vN/Parameters.NFFT*d;
        ExactProfile=LayerInfo(indLayer).n2^2+(LayerInfo(indLayer).n1^2-LayerInfo(indLayer).n2^2)*f_cosrialz(xvet,x1,d1,alpha1);
        [ETEH]=f_EvalProfileFFTH(ExactProfile,vN,Parameters.NHarmonicsTE,d);
    end
    LayerInfo(indLayer).ETEH=ETEH;
    %-- Lossy case
    [kzTEH,kTG_TE2H,VTEH,ITEH]=f_EvalPSWWM_RCWA_TE(k0,LayerInfo(indLayer).NModes,-csiHarmTE,-ky,ETEH);
end
%-- Matrices are already orthogonal; making them orthonormal
Orth=ones(Parameters.NHarmonicsTE,1)*sqrt(diag(VTEH.'*VTE)).';
ITE=ITE./Orth;
ITEH=ITEH./Orth;
VTE=VTE./Orth;
VTEH=VTEH./Orth;
    
LayerInfo(indLayer).kzTE=kzTE;
LayerInfo(indLayer).kTTE=f_psqrt(kTG_TE2);
LayerInfo(indLayer).VTE=VTE.';
LayerInfo(indLayer).ITE=ITE.';
LayerInfo(indLayer).VTEH=VTEH.';
LayerInfo(indLayer).ITEH=ITEH.';
LayerInfo(indLayer).ZG_TE=omega*mu0*kzTE./kTG_TE2;

deltaKronTE=(VTEH/sqrt(d)).'*VTE/sqrt(d);
norm_TE=max(abs(diag(deltaKronTE)*d-ones(LayerInfo(indLayer).NModes,1)));
offdiag=deltaKronTE-diag(diag(deltaKronTE));
ortoP_TE=max(max(abs(offdiag)));

if norm_TE >= 1.e-13
    disp('The diagonal of the eigenvector matrix has not norm 1')
    norm_TE
end
if ortoP_TE >= 1.e-10
    disp('The eigenvector matrix is not orthonormal!')
    ortoP_TE
end

if Options.RebuildProfile==1

    FSMat=exp(-j*x.'*(csiHarmTM-csiHarmTM(1))); %-- Fourier (NOT Floquet) expansion
    figure(13),clf
    axis on
    grid on
    hold on
    RebuiltProfile=FSMat*ETM(:,1);
    plot(x,real(RebuiltProfile),'k')
    ExactProfile=LayerInfo(indLayer).n2^2+(LayerInfo(indLayer).n1^2-LayerInfo(indLayer).n2^2)*f_cosrialz(x,x1,d1,alpha1); 
    plot(x,ExactProfile)

    figure(14),clf
    axis on
    grid on
    hold on
    RebuiltInverseProfile=FSMat*FTM(:,1);
    plot(x,real(RebuiltInverseProfile),'k')
    plot(x,1./ExactProfile)
    pause

end

%-- Determination of PSWW TM modes
[kzTM,kTG_TM2,VTM,ITM]=f_EvalPSWWM_RCWA_TM(k0,LayerInfo(indLayer).NModes,csiHarmTM,ky,ETM,FTM);

if imag(Geometry.n_in)==0 && imag(LayerInfo(indLayer).n1)==0 && imag(LayerInfo(indLayer).n2)==0
    %-- Lossless case
    VTMH=-conj(VTM); %-- "-" sign because of (j*xi*...) instead of (-j*xi*...)
    ITMH=conj(ITM);
else
    %-- Lossy case
    if LayerInfo(indLayer).RollOff==0
        [ETMH,FTMH]=f_EvalRefractiveIndexExpansions(Parameters.NHarmonicsTM,-csiHarmTM,d,alpha1,a1,x1,LayerInfo(indLayer).n1,LayerInfo(indLayer).n2);
    else
        N=Parameters.NFFT/2;
        vN=(-N:N-1); %-- Best choice, to maintain the matlab fft() order!xvet=vN*d/NFFT; %-- Sampling points in natural (spatial) domain
        xvet=vN/Parameters.NFFT*d;
        ExactProfile=LayerInfo(indLayer).n2^2+(LayerInfo(indLayer).n1^2-LayerInfo(indLayer).n2^2)*f_cosrialz(xvet,x1,d1,alpha1);
        [ETMH]=f_EvalProfileFFTH(ExactProfile,vN,Parameters.NHarmonicsTM,d);
        [FTMH]=f_EvalProfileFFTH(1./ExactProfile,vN,Parameters.NHarmonicsTM,d);
    end
    LayerInfo(indLayer).ETMH=ETMH;
    LayerInfo(indLayer).FTMH=FTMH;
    [kzTMH,kTG_TM2H,VTMH,ITMH]=f_EvalPSWWM_RCWA_TM(k0,LayerInfo(indLayer).NModes,-csiHarmTM,-ky,ETMH,FTMH);
end
%-- Matrices are already orthogonal; making them orthonormal
Orth=ones(Parameters.NHarmonicsTM,1)*sqrt(diag(ITMH.'*FTM*ITM)).';
ITM=ITM./Orth;
ITMH=ITMH./Orth;
VTM=VTM./Orth;
VTMH=VTMH./Orth;

LayerInfo(indLayer).kzTM=kzTM;
LayerInfo(indLayer).kTTM=f_psqrt(kTG_TM2);
LayerInfo(indLayer).VTM=VTM.';
LayerInfo(indLayer).ITM=ITM.';
LayerInfo(indLayer).VTMH=VTMH.';
LayerInfo(indLayer).ITMH=ITMH.';
LayerInfo(indLayer).ZG_TM=kTG_TM2./kzTM/(omega*eps0);

%-- For the TM problem, the inner product should be weighted with 1/n^2(x)
deltaKronTM=ITMH.'*FTM*ITM;
normTM=max(abs(diag(deltaKronTM)-ones(LayerInfo(indLayer).NModes,1)));
offdiag=deltaKronTM-diag(diag(deltaKronTM));
ortoPTM=max(max(abs(offdiag)));

if normTM >= 1.e-13
    normTM
end
if ortoPTM >= 1.e-10
    ortoPTM
end

return