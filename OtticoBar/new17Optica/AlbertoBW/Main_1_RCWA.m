
clear
close all
clc





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Geometry of the structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-- Global structure parameters
Geometry.n_in=1; %-- Refractive index left half-space
Geometry.th_in=0; %-- Thickness (z axis) of the left half-space
Geometry.d=1450; %-- Period of the structure
Geometry.th=[240 710   300 840   200 400   950 ]; %-- Thickness of the tooth (z axis)
Geometry.DC=zeros(size(Geometry.th)); %-- Duty cycle of the bar
Geometry.DC([6,7])=0.5;
Geometry.Displacement=zeros(size(Geometry.th)); %-- Displacement of the grating
Geometry.RollOff=zeros(size(Geometry.th)); %-- Roll-off of the raised-cosine profile
Geometry.th_out=0; %-- Thickness (z axis) of the right half-space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Options of the program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Options.RebuildProfile=0; %-- If 1, the layer "n" profile is re-built
Options.PlotGeometry=0; %-- If 1, the geometry is plot
Options.PlotFieldJunctions=0; %-- If 1, junction fields are plotted
Options.PlotTotalField=0; %-- If 1, 3-D field plots are produced
Options.PlotForcedHarmonicContent=0; %-- If 1, 3-D field plots are produced
Options.PlotModeHarmonicContent=0; %-- If 1, 3-D field plots are produced
Options.VerifyPower=0; %-- If 1, the power conservation is verified
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters of the program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Parameters.theta=0;
Parameters.phi=0;
Parameters.vLam=linspace(3100,3500,101);
Parameters.NFFT=2^11;
%-- ALWAYS ODD numbers (generation of Floq. modes)
Parameters.NModes=21; %-- inner hemispace
Parameters.NHarmonicsTE=2*Parameters.NModes+1; %-- Harmonics used to build TE PSWW modes
Parameters.NHarmonicsTM=2*Parameters.NModes+1; %-- Harmonics used to build TM PSWW modes

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Programma
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
theta=Parameters.theta*pi/180;
phi=Parameters.phi*pi/180;

NLam=length(Parameters.vLam);

for indLam=1:NLam
    
    nSi=f_EvalRefractiveIndex_Dispersion(3300/1000,'Si');
    nSiO2=f_EvalRefractiveIndex_Dispersion(3300/1000,'SiO2');
    nPy=f_EvalRefractiveIndex_Dispersion(3300/1000,'Pyrex');
%     nSi=f_EvalRefractiveIndex_Dispersion(Parameters.vLam(indLam)/1000,'Si');
%     nSiO2=f_EvalRefractiveIndex_Dispersion(Parameters.vLam(indLam)/1000,'SiO2');
%     nPy=f_EvalRefractiveIndex_Dispersion(Parameters.vLam(indLam)/1000,'Pyrex');

    Geometry.n1=[nSi nSiO2 nSi nSiO2 nSi nSiO2 nSi ];
    Geometry.n2=[nSi nSiO2 nSi nSiO2 nSi 1     1   ];
    Geometry.n_out=nPy; %-- Refractive index right (outer) half-space

    clear JunctionInfo LayerInfo HalfSpaceInfo
    
    lambda=Parameters.vLam(indLam);
    k0=2*pi/lambda;
    n_in=Geometry.n_in;

    ky=k0*n_in*sin(theta)*sin(phi);
    kx0=k0*n_in*sin(theta)*cos(phi); %-- kx for Floq. harmonic 0 = KBd/d

    [S11,S21,S12,S22,JunctionInfo,LayerInfo,HalfSpaceInfo]=f_EvalSMatrixRCWA(lambda,Parameters,Geometry,Options,kx0,ky);
    
    if Options.PlotFieldJunctions==1
        x=linspace(-Geometry.d/2,Geometry.d/2,1001);
        f_VerifyFieldContinuity(x,k0,lambda,ky,Parameters,JunctionInfo,LayerInfo,HalfSpaceInfo,Geometry)
        
    end
    
    if Options.VerifyPower==1
        indModesLeft=find(imag(HalfSpaceInfo(1).kz)==0);
        indModesRight=find(imag(HalfSpaceInfo(2).kz)==0);
        indModesLeft=[indModesLeft;indModesLeft+Parameters.NModes];
        indModesRight=[indModesRight;indModesRight+Parameters.NModes];
        
        %-- Finding S matrix above cutoff
        S11ac=S11(indModesLeft,indModesLeft);
        S21ac=S21(indModesRight,indModesLeft);
        S12ac=S12(indModesLeft,indModesRight);
        S22ac=S22(indModesRight,indModesRight);
        Sac=[S11ac,S12ac;S21ac,S22ac];
        PowerConvError(indLam)=max(max(abs(Sac*Sac'-eye(size(Sac)))));
    end
    
    if Options.PlotTotalField==1
        z=linspace(0,3500,1001);
        x=linspace(-Geometry.d/2,Geometry.d/2,501)+max(Geometry.Displacement)*2;
%         indMode=1; %-- TE incidence
        indMode=Parameters.NModes(1)+1; %-- TM incidence
        
        a1inc=zeros(2*Parameters.NModes,1);
        a1inc(indMode)=1;
        [Ex,Ey,Ez,Hx,Hy,Hz]=f_PlotField(x,z,lambda,ky,a1inc,Parameters,Geometry,JunctionInfo,LayerInfo,HalfSpaceInfo);
        
        [X,Z]=meshgrid(x,z);
        
        figure(1001),surf(Z,X,abs(Ex),'edgecolor','none'),title('Ex')
        figure(1002),surf(Z,X,abs(Ey),'edgecolor','none'),title('Ey')
        figure(1003),surf(Z,X,abs(Ez),'edgecolor','none'),title('Ez')
        figure(1004),surf(Z,X,abs(Hx),'edgecolor','none'),title('Hx')
        figure(1005),surf(Z,X,abs(Hy),'edgecolor','none'),title('Hy')
        figure(1006),surf(Z,X,abs(Hz),'edgecolor','none'),title('Hz')
        
    end
        
    if Options.PlotForcedHarmonicContent==1
        z=linspace(0,3500,1001);
        x=linspace(-Geometry.d/2,Geometry.d/2,101)+max(Geometry.Displacement)*2;
        
        indMode=1; %-- TE incidence
%         indMode=Parameters.NModes_in(1)+1; %-- TM incidence
        
        a1inc=zeros(2*Parameters.NModes,1);
        a1inc(indMode)=1;
        [VContent,IContent]=f_PlotHarmonicContent_Forced(z,a1inc,Geometry,JunctionInfo,LayerInfo);       
        
        figure
        axis on
        grid on
        hold on
        plot(z,(abs(VContent(:,[indMode]))).^2)
        xlabel('z (nm)')
        ylabel('|V(z)|')
        title(['Tensioni armoniche forzate'])
        
        figure
        axis on
        grid on
        hold on
        plot(z,(abs(IContent(:,[indMode]))).^2)
        xlabel('z (nm)')
        ylabel('|V(z)|')
        title(['Correnti armoniche forzate'])
        
    end
    
    if Options.PlotModeHarmonicContent==1
        z=linspace(0,3200,5001);
        x=linspace(-Geometry.d/2,Geometry.d/2,101)+max(Geometry.Displacement)*2;
        
        indMode=1; %-- TE incidence
%         indMode=Parameters.NModes_in(1)+1; %-- TM incidence

        [VContent]=f_PlotHarmonicContent(z,indMode,Parameters,Geometry,JunctionInfo,LayerInfo);
                
        figure
        axis on
        grid on
        hold on
        plot(z,(abs(VContent(:,[indMode,indMode+3,indMode+1]))).^2./max(abs(VContent(:,[indMode]))).^2*5.31)
        xlabel('z (nm)')
        ylabel('|V(z)|')
        title(['Tensioni armoniche'])
        
    end
    
    S11Tot(:,:,indLam)=S11;
    S12Tot(:,:,indLam)=S12;
    S21Tot(:,:,indLam)=S21;
    S22Tot(:,:,indLam)=S22;
    
    disp(['--- Wavelength ',num2str(indLam),' of ',num2str(length(Parameters.vLam)),' ---'])

end

if Options.VerifyPower==1
    figure(1000001)
    axis on
    grid on
    hold on
    plot(Parameters.vLam,log10(abs(PowerConvError)))
    xlabel('\lambda (nm)')
    ylabel('log10 power conservation error')
end

S11TETERCWA=squeeze(S11Tot(1,1,:));
S21TETERCWA=squeeze(S21Tot(1,1,:));
S11TMTMRCWA=squeeze(S11Tot(Parameters.NModes+1,Parameters.NModes+1,:));
S21TMTMRCWA=squeeze(S21Tot(Parameters.NModes+1,Parameters.NModes+1,:));
S12TETERCWA=squeeze(S12Tot(1,1,:));
S22TETERCWA=squeeze(S22Tot(1,1,:));
S12TMTMRCWA=squeeze(S12Tot(Parameters.NModes+1,Parameters.NModes+1,:));
S22TMTMRCWA=squeeze(S22Tot(Parameters.NModes+1,Parameters.NModes+1,:));

figure(10001)
set(gcf,'position',[560 360 529 588])
subplot(2,1,1)
axis on
grid on
hold on
plot(Parameters.vLam,-log10(1-abs(S11TETERCWA).^2),'linewidth',2)
plot(Parameters.vLam,-log10(1-abs(S11TMTMRCWA).^2),'linewidth',2)
xlim([Parameters.vLam(1),Parameters.vLam(end)])
legend('TE','TM')

return
