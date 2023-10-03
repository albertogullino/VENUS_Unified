
function f_VerifyFieldContinuity(x,k0,lambda,ky,Parameters,JunctionInfo,LayerInfo,HalfSpaceInfo,Geometry)
a1inc=zeros(2*Parameters.NModes,1);
a1inc(1)=1;
% a1inc(16)=1;

for indJunction=1:length(JunctionInfo)
    
%     clear S11L S12L S21L S22L S11C S12C S21C S22C S11R S12R S21R S22R S11CL S12CL S21CL S22CL S11CR S12CR S21CR S22CR
    
    if indJunction==1
        [S11L,S12L,S21L,S22L]=f_CascadeJunctions(JunctionInfo,1,length(JunctionInfo));
        [S11RR,S12RR,S21RR,S22RR]=f_CascadeJunctions(JunctionInfo,2,length(JunctionInfo));
        [S11RL,S12RL,S21RL,S22RL]=f_CascadeJunctions(JunctionInfo,1,1);
        [exL,eyL,ezL,hxL,hyL,hzL]=f_EvalFloquetModes_z(x,k0,HalfSpaceInfo(1).kx,ky,Geometry.d,HalfSpaceInfo(1).n);
        indLayer=JunctionInfo(indJunction).indLayer(2);
        if JunctionInfo(indJunction).TypeLayer(2)==2
            [exR,eyR,ezR,hxR,hyR,hzR]=f_EvalPSWWModes_Floquet(x,lambda,LayerInfo(indLayer),ky,0);
        elseif JunctionInfo(indJunction).TypeLayer(2)==1
            [exR,eyR,ezR,hxR,hyR,hzR]=f_EvalFloquetModes_z(x,k0,HalfSpaceInfo(1).kx,ky,Geometry.d,LayerInfo(indLayer).n);
        end
        kzR=JunctionInfo(1).kzR;
        ZL=diag(JunctionInfo(1).ZinfL);
        ZR=diag(JunctionInfo(1).ZinfR);
        YL=diag(1./JunctionInfo(1).ZinfL);
        YR=diag(1./JunctionInfo(1).ZinfR);
        S11RR=diag(exp(-j.*kzR.*(JunctionInfo(2).zj-JunctionInfo(1).zj)))*S11RR*diag(exp(-j.*kzR.*(JunctionInfo(2).zj-JunctionInfo(1).zj)));
        I=eye(size(S11RR));
        VLinc=sqrt(ZL)*a1inc;
        VLscat=sqrt(ZL)*S11L*a1inc;
        VRscat=sqrt(ZR)*inv(I-S22RL*S11RR)*S21RL*a1inc;
        VRinc=sqrt(ZR)*S11RR*inv(I-S22RL*S11RR)*S21RL*a1inc;
        
    elseif indJunction==length(JunctionInfo)
        [S11L,S12L,S21L,S22L]=f_CascadeJunctions(JunctionInfo,1,length(JunctionInfo)-1);
        [S11C,S12C,S21C,S22C]=f_CascadeJunctions(JunctionInfo,length(JunctionInfo),length(JunctionInfo));
        [S11R,S12R,S21R,S22R]=f_CascadeJunctions(JunctionInfo,1,length(JunctionInfo));
        indLayer=JunctionInfo(indJunction).indLayer(1);
        if LayerInfo(indLayer).flagType==2
            [exL,eyL,ezL,hxL,hyL,hzL]=f_EvalPSWWModes_Floquet(x,lambda,LayerInfo(indLayer),ky,0);
        elseif LayerInfo(indLayer).flagType==1
            [exL,eyL,ezL,hxL,hyL,hzL]=f_EvalFloquetModes_z(x,k0,HalfSpaceInfo(1).kx,ky,Geometry.d,LayerInfo(indLayer).n);
        end    
        [exR,eyR,ezR,hxR,hyR,hzR]=f_EvalFloquetModes_z(x,k0,HalfSpaceInfo(1).kx,ky,Geometry.d,HalfSpaceInfo(2).n);
        
        kzL=JunctionInfo(indJunction).kzL;
        ZL=diag(JunctionInfo(indJunction).ZinfL);
        ZR=diag(JunctionInfo(indJunction).ZinfR);
        YL=diag(1./JunctionInfo(indJunction).ZinfL);
        YR=diag(1./JunctionInfo(indJunction).ZinfR);
        
        S21L=diag(exp(-j.*kzL.*(JunctionInfo(indJunction).zj-JunctionInfo(indJunction-1).zj)))*S21L;
        S22L=diag(exp(-j.*kzL.*(JunctionInfo(indJunction).zj-JunctionInfo(indJunction-1).zj)))*S22L*diag(exp(-j.*kzL.*(JunctionInfo(indJunction).zj-JunctionInfo(indJunction-1).zj)));
        I=eye(size(S11C));
        M=inv(I-S22L*S11C)*S21L;
        VLinc=sqrt(ZL)*M*a1inc;
        VLscat=sqrt(ZL)*S11C*M*a1inc;
        VRscat=sqrt(ZR)*S21R*a1inc;
        VRinc=zeros(size(VRscat));
        
    else
        
        [S11L,S12L,S21L,S22L]=f_CascadeJunctions(JunctionInfo,1,indJunction-1);
        [S11CL,S12CL,S21CL,S22CL]=f_CascadeJunctions(JunctionInfo,indJunction,length(JunctionInfo));
        [S11CR,S12CR,S21CR,S22CR]=f_CascadeJunctions(JunctionInfo,1,indJunction);
        [S11R,S12R,S21R,S22R]=f_CascadeJunctions(JunctionInfo,indJunction+1,length(JunctionInfo));
        indLayerL=JunctionInfo(indJunction).indLayer(1);
        indLayerR=JunctionInfo(indJunction).indLayer(2);
        
        if JunctionInfo(indJunction).TypeLayer(1)==2
            [exL,eyL,ezL,hxL,hyL,hzL]=f_EvalPSWWModes_Floquet(x,lambda,LayerInfo(indLayerL),ky,0);
        elseif JunctionInfo(indJunction).TypeLayer(1)==1
            [exL,eyL,ezL,hxL,hyL,hzL]=f_EvalFloquetModes_z(x,k0,HalfSpaceInfo(1).kx,ky,Geometry.d,LayerInfo(indLayerL).n);
        end        
        
        if JunctionInfo(indJunction).TypeLayer(2)==2
            [exR,eyR,ezR,hxR,hyR,hzR]=f_EvalPSWWModes_Floquet(x,lambda,LayerInfo(indLayerR),ky,0);
        elseif JunctionInfo(indJunction).TypeLayer(2)==1
            [exR,eyR,ezR,hxR,hyR,hzR]=f_EvalFloquetModes_z(x,k0,HalfSpaceInfo(1).kx,ky,Geometry.d,LayerInfo(indLayerR).n);
        end        
        
        kzL=JunctionInfo(indJunction).kzL;
        kzR=JunctionInfo(indJunction).kzR;
        ZL=diag(JunctionInfo(indJunction).ZinfL);
        ZR=diag(JunctionInfo(indJunction).ZinfR);
        YL=diag(1./JunctionInfo(indJunction).ZinfL);
        YR=diag(1./JunctionInfo(indJunction).ZinfR);
        
        S21L=diag(exp(-j.*kzL.*(JunctionInfo(indJunction).zj-JunctionInfo(indJunction-1).zj)))*S21L;
        S22L=diag(exp(-j.*kzL.*(JunctionInfo(indJunction).zj-JunctionInfo(indJunction-1).zj)))*S22L*diag(exp(-j.*kzL.*(JunctionInfo(indJunction).zj-JunctionInfo(indJunction-1).zj)));
        I=eye(size(S22L));
        M=inv(I-S22L*S11CL)*S21L;
        VLinc=sqrt(ZL)*M*a1inc;
        VLscat=sqrt(ZL)*S11CL*M*a1inc;
        S11R=diag(exp(-j.*kzR.*(JunctionInfo(indJunction+1).zj-JunctionInfo(indJunction).zj)))*S11R*diag(exp(-j.*kzR.*(JunctionInfo(indJunction+1).zj-JunctionInfo(indJunction).zj)));
        I=eye(size(S11R));
        M=inv(I-S22CR*S11R)*S21CR;
        VRscat=sqrt(ZR)*M*a1inc;
        VRinc=sqrt(ZR)*S11R*M*a1inc;
    end
   
    VL=VLinc+VLscat;
    IL=YL*(VLinc-VLscat);
    VR=VRinc+VRscat;
    IR=YR*(VRscat-VRinc);
    
    ExL=VL.'*exL;
    ExR=VR.'*exR;
    EyL=VL.'*eyL;
    EyR=VR.'*eyR;
    HxL=IL.'*hxL;
    HxR=IR.'*hxR;
    HyL=IL.'*hyL;
    HyR=IR.'*hyR;
    
    figure(indJunction*10+1),
    subplot(2,1,1)
    axis on
    grid on
    hold on
%     plot(x,abs(ExL),x,abs(ExR),'r')
    plot(x,real(ExL),x,real(ExR),'r')
    title(['Junction ',num2str(indJunction),', Ex'])
    subplot(2,1,2)
    axis on
    grid on
    hold on
%     plot(x,abs(EyL),x,abs(EyR),'r')
    plot(x,real(EyL),x,real(EyR),'r')
    title(['Junction ',num2str(indJunction),', Ey'])
    
    figure(indJunction*10+2),
    subplot(2,1,1)
    axis on
    grid on
    hold on
%     plot(x,abs(HxL),x,abs(HxR),'r')
    plot(x,real(HxL),x,real(HxR),'r')
    title(['Junction ',num2str(indJunction),', Hx'])
    subplot(2,1,2)
    axis on
    grid on
    hold on
%     plot(x,abs(HyL),x,abs(HyR),'r')
    plot(x,real(HyL),x,real(HyR),'r')
    title(['Junction ',num2str(indJunction),', Hy'])

end

return