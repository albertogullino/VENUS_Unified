
function mesh=loadmacro(geom,mesh,mode,T)
%===============================================================================================100
% command preprocessing
if(not(isfield(geom,'mobility_dop'))), geom.mobility_dop='none'; end
if(not(isfield(geom,'ethjflg'))), geom.ethjflg=0; end
if(not(isfield(geom,'dop_t'))), geom.dop_t=[]; end
%===============================================================================================100
% Loading mesh parameters
nt=mesh.nt; % number of triangles
nn=mesh.nn; % number of nodes
%
in1=mesh.triangle(1,:); in2=mesh.triangle(2,:); in3=mesh.triangle(3,:);
it=ismember(mesh.triangle(4,:),find(geom.semiconductor)); % triangles in semiconductor regions
iq=false(1,nn); iq([in1(it) in2(it) in3(it)])=1; mesh.iq=iq;
%
% Loading physical constants
s_LoadConstants
%
% Temperature (K)
Vt=kB.*T./qel;
Ttr=pdeintrp(mesh.node,mesh.triangle(1:4,:),T.'); % T on triangles
%
%==============================================================================================100
% get material properties
mesh.epsxx=zeros(1,nt);
mesh.epsxx_n=zeros(1,nn);
mesh.meffn=zeros(1,nn);
mesh.meffp=zeros(1,nn);
mesh.mobnint_n=zeros(1,nn);
mesh.mobpint_n=zeros(1,nn); % intrinsic mobility
mesh.mobnint_t=zeros(1,nt);
mesh.mobpint_t=zeros(1,nt); % intrinsic mobility
mesh.mobn0_t=zeros(1,nt);
mesh.mobp0_t=zeros(1,nt); % low-field mobility
mesh.mobn0_n=zeros(1,nn);
mesh.mobp0_n=zeros(1,nn); % low-field mobility
mesh.Nc=zeros(1,nn); % DoS
mesh.Nv=zeros(1,nn); % DoS
mesh.affinity=zeros(1,nn); % affinity, eV
mesh.Eg=zeros(1,nn); % energy gap, eV
mesh.nint=zeros(1,nn); % intrinsic concentration, 1/cm^3
mesh.phi_r=zeros(1,nn); % reference potential, V
mesh.DeltaEa=zeros(1,nn); % acceptor impurities activation energy, eV
mesh.DeltaEd=zeros(1,nn); % donor impurities activation energy, eV
mesh.tauc=zeros(1,nn); % capture time for quantum-corrected model, s
mesh.taue=zeros(1,nn); % e time for quantum-corrected model, s
mesh.xmol=zeros(1,nn); % molar fraction in each node
%
for sd=1:geom.nd; % loop on subdomains @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@100
    material=char(geom.material(sd)); % material label
    gvet=geom.gvet{sd};
    it=(mesh.triangle(4,:)==sd); ii=false(1,nn); ii([in1(it) in2(it) in3(it)])=1;
    xtr=mesh.point(1,it); ytr=mesh.point(2,it);
    xnod=mesh.node(1,ii); ynod=mesh.node(2,ii);
    %
    xmol_node=f_EvalMolar(xnod,ynod,gvet); % nodes
    xmol_triangle=f_EvalMolar(xtr,ytr,gvet); % triangles
    %
%    'qui mac', keyboard
    [macro_node]=mw_phmat(material,xmol_node,T(ii),mode);
    [macro_triangle]=mw_phmat(material,xmol_triangle,Ttr(it),mode);
    %
    mesh.epsxx(it)=eps0.*macro_triangle.epsrxx;
    %
    if(geom.semiconductor(sd)), % semiconductor materials
        %
        switch(geom.mobility_dop)
            case('none')
                mesh.mobnint_t(it)=macro_triangle.mobnint; 
                mesh.mobpint_t(it)=macro_triangle.mobpint;
            case('arora_gaas')
                dop=abs(mesh.dop(in1(it))+mesh.dop(in2(it))+mesh.dop(in3(it)))/3;
                [mesh.mobn0(it),mesh.mobp0(it)]=feval(@arora_gaas,T,dop);
                fprintf('Using doping-dependent mobility model arora_gaas for region %d \n',sd)
            otherwise, error('which doping-dependent mobility model?'), end
            %
            mesh.epsxx_n(ii)=eps0.*macro_node.epsrxx;
            mesh.meffn(ii)=macro_node.meffn;
            mesh.meffp(ii)=macro_node.meffp;
            mesh.mobnint_n(ii)=macro_node.mobnint;
            mesh.mobpint_n(ii)=macro_node.mobpint;
            mesh.DeltaEa(ii)=macro_node.DeltaEa;
            mesh.DeltaEd(ii)=macro_node.DeltaEd;
            mesh.Etrap(ii)=macro_node.Etrap;
            mesh.taun(ii)=macro_node.taun;
            mesh.taup(ii)=macro_node.taup;
            mesh.taunQW(ii)=macro_node.taunQW;
            mesh.taupQW(ii)=macro_node.taupQW;            
            mesh.brad(ii)=macro_node.brad;
            mesh.Cnnp(ii)=macro_node.Cnnp;
            mesh.Cppn(ii)=macro_node.Cppn;
            mesh.tauscatn(ii)=macro_node.tauscatn;
            mesh.tauscatn(ii)=macro_node.tauscatn;
            mesh.Nc(ii)=macro_node.Nc;
            mesh.Nv(ii)=macro_node.Nv;
            mesh.affinity(ii)=macro_node.affinity;
            mesh.Eg(ii)=macro_node.Eg;
            mesh.xmol(ii)=xmol_node;
            % mesh.vsurfn(ii)=mesh.arichn(ii).*T^2./(qel*mesh.Nc(ii));
            % mesh.vsurfp(ii)=mesh.arichp(ii).*T^2./(qel*mesh.Nv(ii));
    end
end %@@@@@@@@@@@@@@@@@@@@@@@@@@@@100
%
if(isfield(mode,'BGN') && mode.BGN == 1)
%     IBTJ=unique(cell2mat(mesh.IBTJ));
%     mesh.Eg(IBTJ)= mesh.Eg(IBTJ) - 3.5e-8.*mesh.dop_d(IBTJ).^(1/3) - 3e-8.*mesh.dop_a(IBTJ).^(1/3); % eV
    if isfield(mode,'X_D')
        mesh.Eg= mesh.Eg - mode.X_D*1e-8.*(mesh.dop_d*mode.CarrierNorm).^(1/3) - mode.X_A*1e-8.*(mesh.dop_a*mode.CarrierNorm).^(1/3); % eV
    else
        mesh.Eg= mesh.Eg - 3.5e-8.*(mesh.dop_d*mode.CarrierNorm).^(1/3) - 3e-8.*(mesh.dop_a*mode.CarrierNorm).^(1/3); % eV
    end
end
%
% if mode.flgBTJ==1  % flag to TJ presence
%     IBTJ=unique(cell2mat(mesh.IBTJ));
%     mesh.DeltaEa(IBTJ)=mode.EaTJ;
%     mesh.DeltaEd(IBTJ)=mode.EdTJ;
% end
% reference potential, V
mesh.phi_r(iq)=-mesh.affinity(iq)-1/2.*mesh.Eg(iq)+1/2.*Vt(iq).*log(mesh.Nv(iq)./mesh.Nc(iq));
iq_r=find(iq); iq_r=iq_r(1);
mesh.phi_rr=mesh.phi_r(iq_r);
mesh.phi_r(iq)=mesh.phi_r(iq)-mesh.phi_r(iq_r);
%
for ic=1:geom.nc; ii=(mesh.contact==ic); jj=find(ii&iq);
    if(length(jj)>0), mesh.phi_r(ii)=mesh.phi_r(jj(1)); end, end
%
% intrinsic concentration, 1/cm^3
mesh.nint(iq)=sqrt(mesh.Nc(iq).*mesh.Nv(iq)).*exp(-mesh.Eg(iq)./(2*Vt(iq)));
mesh.T=T;
%
% corrections of mobility by doping
[mesh] = f_EvalDopingMobility(mesh,mode);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Temperature-dependent quantum-correction related quantities
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(isfield(geom,'QWorientation'))
    CarrierNorm2D=mode.CarrierNorm2D;
    for indQW=1:mesh.NMQW
        QWcells=mesh.MQWcell{indQW};
        TQW=T(mesh.inMQW{indQW});
        % [macroQW]=mw_phmat(geom.material{QWcells},geom.QWxmol(indQW),TQW);
        xmol=imag(geom.QWxmol{indQW});
        [macroQW]=mw_phmat(geom.material{QWcells},xmol,TQW,mode);
        mesh.Nc2D{indQW}=4*pi*mesh.meffnMQW{indQW}*m0*kB*TQW/h^2/10000/CarrierNorm2D;
        mesh.Nv2D{indQW}=4*pi*mesh.meffpMQW{indQW}*m0*kB*TQW/h^2/10000/CarrierNorm2D;
        if(isfield(macroQW,'mobnint'))
            mesh.mobnQW{indQW}=0.5*(macroQW.mobnint(1:end-1)+macroQW.mobnint(2:end));
            mesh.mobpQW{indQW}=0.5*(macroQW.mobpint(1:end-1)+macroQW.mobpint(2:end));
            mesh.tauscatnMQW{indQW}=macroQW.tauscatn;
            mesh.tauscatpMQW{indQW}=macroQW.tauscatp;
        end
        if abs(mesh.xmol_barrierLUT-mesh.xmol(mesh.inMQW{1}(1)))>0
            mesh.DeltaEcQW{indQW}=mesh.DeltaEgLUT*mesh.xmol(mesh.inMQW{indQW}(1))*mesh.QcLUT;
            mesh.DeltaEvQW{indQW}=mesh.DeltaEgLUT*mesh.xmol(mesh.inMQW{indQW}(1))*(1-mesh.QcLUT);
        end
    end
end