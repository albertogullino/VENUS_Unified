function [ex,ey,ez,hx,hy,hz]=f_EvalPSWWModes_Floquet(x,lambda,LayerInfo,ky,flagAdj)
%-- questa funzione calcola i modi blabla nella base di Floquet. In
%pratica, ricostruisce usando gli autovettori salvati nella struttura
%LayerInfo.

f_LoadConstants

%-- Number of Floquet modes to represent the grating modes
Nmodi=size(LayerInfo.VTE,1);
d=LayerInfo.d1+LayerInfo.d2;
k0=2*pi/lambda;
omega=k0*Clight;

if flagAdj==0
    
    csiHarmTE=LayerInfo.csiHarmTE;
    csiHarmTM=LayerInfo.csiHarmTM;
    
    FloqMatTE=exp(-j*csiHarmTE.'*x)./sqrt(d); %-- Ricostruzione campi/tutto
    FloqMatTM=exp(-j*csiHarmTM.'*x)./sqrt(d); %-- Ricostruzione campi/tutto
    
    %-- formule pagina 10 report Orta (TE Modes)
    
    eyTE=(LayerInfo.VTE*FloqMatTE);
    exTE=zeros(size(eyTE));
    ezTE=-((ky./LayerInfo.kzTE)*ones(1,length(x))).*(LayerInfo.VTE*FloqMatTE);
    hxTE=-(LayerInfo.VTE*FloqMatTE);
    hyTE=((omega.*mu0.*ky./(LayerInfo.kTTE).^2)*ones(1,length(x))).*(LayerInfo.ITE*FloqMatTE);
    hzTE=(LayerInfo.ZG_TE*ones(1,length(x))).*(LayerInfo.ITE*FloqMatTE);
    
    %-- TM modes
    d1=LayerInfo.d1;
    n1=LayerInfo.n1;
    n2=LayerInfo.n2;
    alpha1=LayerInfo.RollOff;
    x1=LayerInfo.Displacement;
    [n]=f_EvalRefractiveIndex(x,d1,n1,n2,alpha1,x1);
%     figure,plot(x,n.^2)
%     return
    VN2=ones(Nmodi,1)*(n.^2);
    
    exTM=(LayerInfo.ITM*FloqMatTM)./VN2;
    %-- NOTA CHE il ./VN2 è per via della definizione, (2.4), di V'
    eyTM=-((ky./(LayerInfo.kTTM.^2./(omega.*eps0)))*ones(1,length(x))).*(LayerInfo.VTM*FloqMatTM)./VN2;
    ezTM=-((1./(LayerInfo.ZG_TM))*ones(1,length(x))).*(LayerInfo.VTM*FloqMatTM)./VN2;
    hxTM=zeros(size(exTM));
    hyTM=(LayerInfo.ITM*FloqMatTM);
    hzTM=-((ky./(LayerInfo.kzTM))*ones(1,length(x))).*(LayerInfo.ITM*FloqMatTM);
    
elseif flagAdj==1
    
    csiHarmTE=-LayerInfo.csiHarmTE;
    csiHarmTM=-LayerInfo.csiHarmTM;
    ky = -ky;
    
    FloqMatTE=exp(-j*csiHarmTE.'*x)./sqrt(d); %-- Ricostruzione campi/tutto
    FloqMatTM=exp(-j*csiHarmTM.'*x)./sqrt(d); %-- Ricostruzione campi/tutto
    
    %-- formule pagina 10 report Orta (TE Modes)
    
    eyTE=(LayerInfo.VTEH*FloqMatTE);
    exTE=zeros(size(eyTE));
    ezTE=-(((ky)./LayerInfo.kzTE)*ones(1,length(x))).*(LayerInfo.VTEH*FloqMatTE);
    hxTE=-(LayerInfo.VTEH*FloqMatTE);
    hyTE=((omega.*mu0.*(-ky)./(LayerInfo.kTTE).^2)*ones(1,length(x))).*(LayerInfo.ITEH*FloqMatTE);
    hzTE=(LayerInfo.ZG_TE*ones(1,length(x))).*(LayerInfo.ITEH*FloqMatTE);
    
    %-- TM modes
    d1=LayerInfo.d1;
    n1=LayerInfo.n1;
    n2=LayerInfo.n2;
    [n]=f_EvalRefractiveIndex(x,d1,n1,n2,alpha1,x1);
    
    VN2=ones(Nmodi,1)*(n.^2);
    
    exTM=(LayerInfo.ITMH*FloqMatTM)./VN2;
    %-- NOTA CHE il ./VN2 è per via della definizione, (2.4), di V'
    eyTM=-(((-ky)./(LayerInfo.kTTM.^2./(omega.*eps0)))*ones(1,length(x))).*(LayerInfo.VTMH*FloqMatTM)./VN2;
    ezTM=-((1./(LayerInfo.ZG_TM))*ones(1,length(x))).*(LayerInfo.VTMH*FloqMatTM)./VN2;
    hxTM=zeros(size(exTM));
    hyTM=(LayerInfo.ITMH*FloqMatTM);
    hzTM=-((ky./(LayerInfo.kzTM))*ones(1,length(x))).*(LayerInfo.ITMH*FloqMatTM);
end

ex=[exTE;exTM];
ey=[eyTE;eyTM];
ez=[ezTE;ezTM];
hx=[hxTE;hxTM];
hy=[hyTE;hyTM];
hz=[hzTE;hzTM];

return
