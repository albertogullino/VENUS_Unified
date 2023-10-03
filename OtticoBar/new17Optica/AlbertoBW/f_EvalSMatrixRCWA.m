
function [S11,S21,S12,S22,JunctionInfo,LayerInfo,HalfSpaceInfo,PowerConvError]=f_EvalSMatrixRCWA(lambda,Parameters,Geometry,Options,kx0,ky)

f_LoadConstants

NLayers=length(Geometry.n1);
k0=2*pi/lambda;
n_in=Geometry.n_in;
n_out=Geometry.n_out;

vp_in_temp = 1:(Parameters.NModes-1)/2;
vp_in(2*vp_in_temp) = -vp_in_temp;
vp_in(2*vp_in_temp+1) = vp_in_temp;

vp_out_temp = 1:(Parameters.NModes-1)/2;
vp_out(2*vp_out_temp) = -vp_out_temp;
vp_out(2*vp_out_temp+1) = vp_out_temp;

%-- Floquet mode parameters go in HalfSpaceInfo structure
%-- HalfSpaceInfo(1): left hemispace
HalfSpaceInfo(1).vp=vp_in;
HalfSpaceInfo(1).kx=kx0 + 2*pi*HalfSpaceInfo(1).vp/Geometry.d;
HalfSpaceInfo(1).kTF=f_psqrt(HalfSpaceInfo(1).kx.^2+ky^2).';
HalfSpaceInfo(1).kz=f_psqrt((k0*n_in)^2 - HalfSpaceInfo(1).kTF.^2);
HalfSpaceInfo(1).YinfTE=HalfSpaceInfo(1).kz./(k0*Clight*mu0);
HalfSpaceInfo(1).YinfTM=k0*n_in^2./(HalfSpaceInfo(1).kz*Clight*mu0);
HalfSpaceInfo(1).n=n_in;
%-- HalfSpaceInfo(2): right hemispace
HalfSpaceInfo(2).vp=vp_out;
HalfSpaceInfo(2).kx=kx0 + 2*pi*HalfSpaceInfo(2).vp/Geometry.d;
HalfSpaceInfo(2).kTF=f_psqrt(HalfSpaceInfo(2).kx.^2+ky^2).';
HalfSpaceInfo(2).kz=f_psqrt((k0*n_out)^2 - HalfSpaceInfo(2).kTF.^2);
HalfSpaceInfo(2).YinfTE=HalfSpaceInfo(2).kz./(k0*Clight*mu0);
HalfSpaceInfo(2).YinfTM=k0*n_out^2./(HalfSpaceInfo(2).kz*Clight*mu0);
HalfSpaceInfo(2).n=n_out;

%-- Loading data in LayerInfo structure
for indLayer=1:NLayers
    if abs(Geometry.n1(indLayer)-Geometry.n2(indLayer))<1e-12
        LayerInfo(indLayer).flagType=1; %-- Type 1 is "homogeneous"
        LayerInfo(indLayer).n=Geometry.n1(indLayer);
        LayerInfo(indLayer).th=Geometry.th(indLayer);
        LayerInfo(indLayer).NModes=Parameters.NModes;
        kzTE=f_psqrt((k0*LayerInfo(indLayer).n).^2-HalfSpaceInfo(1).kTF.^2);
        kzTM=f_psqrt((k0*LayerInfo(indLayer).n).^2-HalfSpaceInfo(1).kTF.^2);
        LayerInfo(indLayer).kzTE=kzTE;
        LayerInfo(indLayer).kzTM=kzTM;
        LayerInfo(indLayer).ZG_TE=(k0*Clight*mu0)./kzTE;
        LayerInfo(indLayer).ZG_TM=kzTM.*Clight.*mu0./(k0.*LayerInfo(indLayer).n.^2);
    else
        %-- LayerInfo structure is for inner (tooth) layers!
        LayerInfo(indLayer).flagType=2; %-- Type 2 is "grating"
        LayerInfo(indLayer).n1=Geometry.n1(indLayer);
        LayerInfo(indLayer).n2=Geometry.n2(indLayer);
        LayerInfo(indLayer).d1=Geometry.DC(indLayer)*Geometry.d;
        LayerInfo(indLayer).d2=(1-Geometry.DC(indLayer))*Geometry.d;
        LayerInfo(indLayer).th=Geometry.th(indLayer);
        LayerInfo(indLayer).NModes=Parameters.NModes;
        LayerInfo(indLayer).RollOff=Geometry.RollOff(indLayer);
        LayerInfo(indLayer).Displacement=Geometry.Displacement(indLayer);
        [LayerInfo]=f_EvalPSWW_Modes(lambda,indLayer,Parameters,Geometry,LayerInfo,Options,kx0,ky);
    end
end

%-- Solving junction problems and saving in JunctionInfo structure

%-- First junction
indJunction=1;
indLayer_right=1;
ZL=[1./HalfSpaceInfo(1).YinfTE;1./HalfSpaceInfo(1).YinfTM];
ZR=[LayerInfo(indLayer_right).ZG_TE;LayerInfo(indLayer_right).ZG_TM];
JunctionInfo(indJunction).zj=Geometry.th_in;
JunctionInfo(indJunction).kzL=[HalfSpaceInfo(1).kz;HalfSpaceInfo(1).kz];
JunctionInfo(indJunction).ZinfL=ZL;
JunctionInfo(indJunction).ZinfR=ZR;
JunctionInfo(indJunction).kzR=[LayerInfo(indLayer_right).kzTE;LayerInfo(indLayer_right).kzTM];
JunctionInfo(indJunction).indLayer=[1,indLayer_right]; %-- Sx e Dx

if LayerInfo(indLayer_right).flagType==1
    [S11j,S21j,S12j,S22j]=f_EvalMLJunction(ZL,ZR);
    JunctionInfo(indJunction).TypeLayer=[1,1];
elseif LayerInfo(indLayer_right).flagType==2
    [S11j,S21j,S12j,S22j]=f_EvalMMTJunction1(LayerInfo,indLayer_right,Parameters,HalfSpaceInfo,ky);
    JunctionInfo(indJunction).TypeLayer=[1,2];
end

JunctionInfo(indJunction).S11=S11j;
JunctionInfo(indJunction).S21=S21j;
JunctionInfo(indJunction).S12=S12j;
JunctionInfo(indJunction).S22=S22j;

if NLayers>1
    %-- Junctions from 2 to NLayers
    for indJunction=2:NLayers
        indLayer_left=indJunction-1;
        indLayer_right=indJunction;
        ZL=[LayerInfo(indLayer_left).ZG_TE;LayerInfo(indLayer_left).ZG_TM];
        ZR=[LayerInfo(indLayer_right).ZG_TE;LayerInfo(indLayer_right).ZG_TM];
        
        JunctionInfo(indJunction).zj=JunctionInfo(indJunction-1).zj+Geometry.th(indLayer_left);
        JunctionInfo(indJunction).ZinfL=ZL;
        JunctionInfo(indJunction).kzL=[LayerInfo(indLayer_left).kzTE;LayerInfo(indLayer_left).kzTM];
        JunctionInfo(indJunction).ZinfR=ZR;
        JunctionInfo(indJunction).kzR=[LayerInfo(indLayer_right).kzTE;LayerInfo(indLayer_right).kzTM];
        JunctionInfo(indJunction).TypeLayer=[LayerInfo(indLayer_left).flagType,LayerInfo(indLayer_right).flagType];
        JunctionInfo(indJunction).indLayer=[indLayer_left,indLayer_right];
        
        if LayerInfo(indLayer_left).flagType==1 && LayerInfo(indLayer_right).flagType==1
            [S11j,S21j,S12j,S22j]=f_EvalMLJunction(ZL,ZR);
        elseif LayerInfo(indLayer_left).flagType==2 && LayerInfo(indLayer_right).flagType==2
            %-- Choosing the most suitable mode-matching formulation
            deltanL=abs(LayerInfo(indLayer_left).n1-LayerInfo(indLayer_left).n2);
            deltanR=abs(LayerInfo(indLayer_right).n1-LayerInfo(indLayer_right).n2);
            if deltanL<deltanR
                [S11j,S21j,S12j,S22j]=f_EvalMMTJunctionLayer_RL(lambda,LayerInfo,indLayer_right,Parameters,ky);
            else
                [S11j,S21j,S12j,S22j]=f_EvalMMTJunctionLayer_LR(lambda,LayerInfo,indLayer_right,Parameters,ky);
            end
        elseif LayerInfo(indLayer_left).flagType==1 && LayerInfo(indLayer_right).flagType==2
            [S11j,S21j,S12j,S22j]=f_EvalMMTJunctionFG(LayerInfo,indLayer_right,Parameters,HalfSpaceInfo,ky);
        elseif LayerInfo(indLayer_left).flagType==2 && LayerInfo(indLayer_right).flagType==1
            [S11j,S21j,S12j,S22j]=f_EvalMMTJunctionGF(LayerInfo,indLayer_right,Parameters,HalfSpaceInfo,ky);
        end
        
        JunctionInfo(indJunction).S11=S11j;
        JunctionInfo(indJunction).S21=S21j;
        JunctionInfo(indJunction).S12=S12j;
        JunctionInfo(indJunction).S22=S22j;
        
    end
end

%-- Ending junction
indJunction=indJunction+1;
indLayer_left=indLayer_right;
JunctionInfo(indJunction).zj=JunctionInfo(indLayer_left).zj+Geometry.th(indLayer_left);
JunctionInfo(indJunction).ZinfL=[LayerInfo(indLayer_left).ZG_TE;LayerInfo(indLayer_left).ZG_TM];
JunctionInfo(indJunction).kzL=[LayerInfo(indLayer_left).kzTE;LayerInfo(indLayer_left).kzTM];
JunctionInfo(indJunction).ZinfR=[1./HalfSpaceInfo(2).YinfTE;1./HalfSpaceInfo(2).YinfTM];
JunctionInfo(indJunction).kzR=[HalfSpaceInfo(2).kz;HalfSpaceInfo(2).kz];
JunctionInfo(indJunction).TypeLayer=[LayerInfo(indLayer_left).flagType,1];

if LayerInfo(indLayer_left).flagType==1
    ZL=JunctionInfo(indJunction).ZinfL;
    ZR=JunctionInfo(indJunction).ZinfR;
    [S11j,S21j,S12j,S22j]=f_EvalMLJunction(ZL,ZR);
elseif LayerInfo(indLayer_left).flagType==2
    [S11j,S21j,S12j,S22j]=f_EvalMMTJunctionEnd(LayerInfo,indLayer_left,Parameters,HalfSpaceInfo,ky);
end

JunctionInfo(indJunction).S11=S11j;
JunctionInfo(indJunction).S21=S21j;
JunctionInfo(indJunction).S12=S12j;
JunctionInfo(indJunction).S22=S22j;
JunctionInfo(indJunction).indLayer=[indLayer_left,2];

if Options.PlotGeometry==1
    f_PlotGeometry(JunctionInfo,LayerInfo,Geometry)
end

[S11,S12,S21,S22]=f_CascadeJunctions(JunctionInfo,1,length(JunctionInfo));

KzIn=[HalfSpaceInfo(1).kz;HalfSpaceInfo(1).kz];
KzOut=[HalfSpaceInfo(2).kz;HalfSpaceInfo(2).kz];

S11=diag(exp(-j.*KzIn.*Geometry.th_in))*S11*diag(exp(-j.*KzIn.*Geometry.th_in));
S12=diag(exp(-j.*KzIn.*Geometry.th_in))*S12*diag(exp(-j.*KzOut.*Geometry.th_out));
S21=diag(exp(-j.*KzOut.*Geometry.th_out))*S21*diag(exp(-j.*KzIn.*Geometry.th_in));
S22=diag(exp(-j.*KzOut.*Geometry.th_out))*S22*diag(exp(-j.*KzOut.*Geometry.th_out));
   
return