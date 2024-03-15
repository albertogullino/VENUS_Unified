
function [Ban,mesh]=f_ComputeQWSubbands(mesh,mode)

%==========================================================================
% Definition of FEM grid and structure parameters
%==========================================================================

vz_offset=0;
% vz_offset=[-mesh.Lz,mesh.Lz]; % 2 QWs case
% vz_offset=[-2*mesh.Lz,0,2*mesh.Lz]; % 3 QWs case
s_LoadConstants
s_GaAsConstants

mesh.ne=mesh.nn-1; % Number of elements
mesh.x=linspace(0,mesh.L,mesh.nn); % Node coordinates, m
mesh.le=mesh.x(2:end)-mesh.x(1:end-1); % Edge length, m
mesh.xc=(mesh.x(1:end-1)+mesh.x(2:end))/2; % Center points, m

xmol_well=mesh.xmol_well;
xmol_barrier=mesh.xmol_barrier;

Delc=mesh.Qc*xmol_barrier;
Delh=(1-mesh.Qc)*xmol_barrier;

C0=Delc*(mesh.DeltaEg); % Conduction band offset, eV
V0=-Delh*(mesh.DeltaEg); % Valence band offset, eV

mesh.C0=C0;
mesh.V0=V0;
mesh.ecb=C0*ones(1,mesh.ne); % Conduction band offset vector
mesh.evb=V0*ones(1,mesh.ne); % Valence band offset vector
mesh.xmol=xmol_barrier*ones(1,mesh.ne);



for indQW=1:length(vz_offset)
    ii=((mesh.xc>=(mesh.L/2-mesh.Lz/2+vz_offset(indQW))) & (mesh.xc<=(mesh.L/2+mesh.Lz/2+vz_offset(indQW))));
    mesh.evb(ii)=0;
    mesh.ecb(ii)=0;
    mesh.xmol(ii)=xmol_well;
end

if mode.iplot==1
    %
    figure(1)
    set(gcf,'Position',[90 524 560 420])
    axis on
    grid on
    hold on
    box on
    plot(mesh.xc*1e9,mesh.ecb+mesh.Eg,'bo')
    plot(mesh.xc*1e9,mesh.evb,'ro')
    xlabel('z (nm)')
    ylabel('Valence band structure, eV')
    drawnow
end

kx=linspace(0,mesh.max_k,mesh.num_kvectors); % kx grid points, angstrom
kx(1)=1e-5; % k grid points along x, angstrom

ky=zeros(1,mesh.num_kvectors); % k grid points along y, angstrom
kgrid=sqrt(kx.^2+ky.^2)*1e10; % k vectors grid, conversion to 1/m

SBC=zeros(mesh.ncinit,mesh.num_kvectors); % CB subbands
XVC=zeros(4*mesh.nn,mesh.ncinit,mesh.num_kvectors); % CB eigenfunctions
SBV=zeros(mesh.nvinit,mesh.num_kvectors); % VB subbands
XVV=zeros(4*mesh.nn,mesh.nvinit,mesh.num_kvectors); % VB eigenfunctions

%==========================================================================
% Computations at k=0
%==========================================================================
mesh.targetV=max(mesh.evb); % target for eigenvalue solver, valence band
mesh.targetC=min(mesh.ecb)+1.5; % target for eigenvalue solver, valence band
[lmbVtmp,xvVtmp,lmbCtmp,xvCtmp]=solve_kp44_AlGaAs(mesh,0,0);

xvC=reshape(xvCtmp,4,mesh.nn,mesh.ncinit);
xvV=reshape(xvVtmp,4,mesh.nn,mesh.nvinit);

% Due to the orthogonality of the orbital basis functions, the
% contributions from the 4x4 k.p results to the dipole matrix element are
% the same as in the 2x2. So, most contributions to the eigenfunctions can
% be ignored.
xvC=squeeze(xvC(1,:,:)); % first index: electrons
xvV1=squeeze(xvV(2,:,:)); % second index: heavy holes
xvV2=squeeze(xvV(3,:,:)); % third index: light holes

for indSB=1:mesh.ncinit
    vnormC(indSB)=norm(xvC(:,indSB));
end

for indSB=1:mesh.nvinit
    vnorm1(indSB)=norm(xvV1(:,indSB));
    vnorm2(indSB)=norm(xvV2(:,indSB));
end

if(vnorm1(1)>1e-5)
    mesh.indhh=find(vnorm1>1e-5);
    mesh.indlh=find(vnorm2>1e-5);
elseif(vnorm2(1)>1e-5)
    mesh.indhh=find(vnorm2>1e-5);
    mesh.indlh=find(vnorm1>1e-5);
end




disp('Computing subbands')

for ik=1:mesh.num_kvectors
    fprintf('k, 1/A %e %e\n',[kx(ik) ky(ik)])
    [lmbV,xvV,lmbC,xvC]=solve_kp44_AlGaAs(mesh,kx(ik),ky(ik));
    
    % Normalization of the valence-like eigenfunction
    for ib=1:mesh.nvinit
        xv=xvV(:,ib);
        [~,ind]=max(abs(xv));
        xv=xv*sign(xv(ind));
        f=reshape(xv,4,mesh.nn);
        g=sum(abs(f).^2,1);
        integ=trapz(mesh.x,g);
        xv=xv/sqrt(integ);
        integ=trapz(mesh.x,g);
        XVV(:,ib,ik)=xv;
        SBV(ib,ik)=lmbV(ib);
    end
    
    % Normalization of the conduction-like eigenfunction
    for ib=1:mesh.ncinit
        xv=xvC(:,ib);
        [~,ind]=max(abs(xv));
        xv=xv*sign(xv(ind));
        f=reshape(xv,4,mesh.nn);
        g=sum(abs(f).^2,1);
        integ=trapz(mesh.x,g);
        xv=xv/sqrt(integ);
        integ=trapz(mesh.x,g);
        XVC(:,ib,ik)=xv;
        SBC(ib,ik)=lmbC(ib);
    end
    
    mesh.targetC=min(lmbC);
    mesh.targetV=max(lmbV);
    
end

xxVC=reshape(XVC,4,mesh.nn,mesh.ncinit,mesh.num_kvectors);
xxVV=reshape(XVV,4,mesh.nn,mesh.nvinit,mesh.num_kvectors);

XVC=squeeze(xxVC(1,:,:,:));

%==========================================================================
% Discard subbands greater than the AlGaAs barrier
%==========================================================================
SBV=-real(SBV);
SBC=real(SBC);
%
indBound_V=find(SBV(:,1)<=-V0);
indBound_C=find(SBC(:,1)<=+(C0+mesh.Eg));
mesh.ncb=length(indBound_C);
mesh.nvb=length(indBound_V);

XVV=XVV(:,indBound_V,:);
SBV=SBV(indBound_V,:);
XVC=XVC(:,indBound_C,:);
SBC=SBC(indBound_C,:);

%==========================================================================
% Sorting relevant (for dipole) eigenfunctions
%==========================================================================
dz=diff(mesh.x(1:2));

XV1=XVV(2:4:end,:,:); % second index: heavy holes
XV2=XVV(3:4:end,:,:); % third index: light holes
%
nk=size(XV1,3);

for indk=2:nk
    
    xv1_old = squeeze(XV1(:,:,indk-1));
    xv1_new = squeeze(XV1(:,:,indk));
    xv2_old = squeeze(XV2(:,:,indk-1));
    xv2_new = squeeze(XV2(:,:,indk));
    
    ovp=xv1_old'*dz*xv1_new+xv2_old'*dz*xv2_new;
    
    [val,indsort]=max(abs(ovp),[],2);
    
    ovp2=ovp(:,indsort);
    XV1(:,:,indk)=XV1(:,indsort,indk);
    XV2(:,:,indk)=XV2(:,indsort,indk);
    XVV(:,:,indk)=XVV(:,indsort,indk);
    SBV(:,indk)=SBV(indsort,indk);
    
end

%======================================================================
% Computing dipole momentum
%======================================================================
Eg=qel*mesh.Eg; % bandgap, J
Delta=qel*mesh.Delta; % spin-orbit coupling, J
mc=mesh.meffn*m0; % conduction band effective mass, kg

Lz=mesh.Lz;
% From 1990 Bava, Debernardi, Lin report, p. 20
R2=(qel.*hbar/2)^2*(Eg.*(Eg+Delta)/(Eg+2/3*Delta))*(1/mc-1/m0);

% Normalization constant for gain and spontaneous emission and
% recombination
ti=4/(eps0*Nb^2) * 1/(2*pi)^2 * (2*pi)/Lz/hbar;

kdk=[0 kgrid(2:end).*diff(kgrid)]';
kdk(end)=kdk(end)/2;

nk=size(XV1,3); % number of points in k grid
nvl=size(XV1,2); % number of light hole bands
nvh=size(XV2,2); % number of heavy hole bands
nc=size(XVC,2); % number of conduction bands

XV1=reshape(XV1,mesh.nn,nk*nvh);
XV2=reshape(XV2,mesh.nn,nk*nvl);
XVC=reshape(XVC,mesh.nn,nk*nc);

% heavy-hole / electron overlap
SHtmp=XV1'*dz*XVC; % this is s1 from 1990Bava report, p. 20
SHtmp=reshape(SHtmp,nvh,nk,nc,nk);
SH=zeros(nvh,nk,nc);
for indnvh=1:nvh
    for indnc=1:nc
        SH(indnvh,:,indnc)=diag(squeeze(SHtmp(indnvh,:,indnc,:)));
    end
end

% light-hole / electron overlap
SLtmp=XV2'*dz*XVC; % this is s2 from 1990Bava report, p. 20
SLtmp=reshape(SLtmp,nvl,nk,nc,nk);
SL=zeros(nvl,nk,nc);
for indnvl=1:nvl
    for indnc=1:nc
        SL(indnvl,:,indnc)=diag(squeeze(SLtmp(indnvl,:,indnc,:)));
    end
end

sih=size(SH);
SHt=[]; SLt=[]; ECV=[];

for ke=1:nc
    SHt=[SHt; SH(:,:,ke).^2];
    SLt=[SLt; SL(:,:,ke).^2];
    ECV=[ECV; SBV+ones(sih(1),1)*SBC(ke,:)];
end

%======================================================================
% Computing potentials for many-body effects
%======================================================================
XVC=reshape(XVC,mesh.nn,mesh.ncb,nk);
XVC=squeeze(XVC(:,:,1));
XV1d=reshape(XV1,mesh.nn,mesh.nvb,nk);
XV2d=reshape(XV2,mesh.nn,mesh.nvb,nk);
XV1=squeeze(XV1d(:,:,1));
XV2=squeeze(XV2d(:,:,1));
XVV=XV1+XV2;
Q=qel;
nbc=mesh.ncb;
nbv=mesh.nvb;
Vcc = zeros(1,length(kgrid));
Vvv = zeros(1,length(kgrid));
[X,Y] = meshgrid(mesh.x,mesh.x); XY = abs(Y-X);
dz = mesh.L/mesh.ne;

% Setting q value to prevent singularities at denominators
dk = kgrid(3) - kgrid(2);
qs = 1e7;

% Gauss-Legendre quadrature nodes for efficient overlap computation
% for potentials
Nk=21;
Nth=31;
kmax = kgrid(ceil(end/2));

% radial k grid (krho)
[nodes,weights]=quadad('legen',1,Nk);
vkrho=kmax/2*nodes+kmax/2;
wkrho=kmax/2*weights;

% angular k grid (theta)
[nodes,weights]=quadad('legen',1,Nth);
vktheta=2*pi/2*nodes+2*pi/2;
wktheta=2*pi/2*weights;

[VKTHETA,VKRHO]=meshgrid(vktheta,vkrho);

for ic = 1:mesh.ncb
    f14 = abs(XVC(:,ic).').^2;
    PPEcV{ic} = spline(kgrid,SBC(ic,:));
    for iq = 1:length(kgrid)
        qD = kgrid(iq) + qs;
        % Unscreened conduction potential, Veprek PhD, (5.2)
        Vcc(iq) = qel^2/(2*eps0*ER_STA*qD)*(dz^2)*(f14*exp(-qD*XY)*f14.');
    end
    % Spline coefficients for conduction potential
    PPVccV{ic} = spline(kgrid,Vcc);
end

for iv = 1:mesh.nvb
    f23 = abs(XVV(:,iv).').^2;
    PPEvV{iv} = spline(kgrid,SBV(iv,:));
    for iq = 1:length(kgrid)
        qD = kgrid(iq) + qs;
        % Unscreened valence potential, Veprek PhD, (5.2)
        Vvv(iq) = qel^2/(2*eps0*ER_STA*qD)*(dz^2)*(f23*exp(-qD*XY)*f23.');
    end
    % Spline coefficients for valence potential
    PPVvvV{iv} = spline(kgrid,Vvv);
end

%======================================================================
% Saving data from band, dipole and potential calculations
%======================================================================
Ban.vktheta=vktheta;
Ban.wkrho=wkrho;
Ban.wktheta=wktheta;
Ban.VKRHO=VKRHO;
Ban.VKTHETA=VKTHETA;
Ban.kdk=kdk;
Ban.kgrid=kgrid;
Ban.PPEcV=PPEcV;
Ban.PPEvV=PPEvV;
Ban.PPVccV=PPVccV;
Ban.PPVvvV=PPVvvV;
Ban.Nb=Nb;
Ban.nvl=nvl;
Ban.nvh=nvh;
Ban.XVC=XVC;
Ban.XVV=XVV;

if(mesh.parabolic==1)
    
    h2m0=hbar^2/(2*m0);

    % For GaAs only !!!!!! -.-
    meffn=0.082; % best fit
    %
    gamma1=6.85;
    gamma2=2.10;
    meffp=(1/(gamma1-2*gamma2))*1.5; % heavy hole mass only
    %
    lmbC=SBC(:,1);
    SBCtmp=lmbC*ones(1,nk)+ones(nc,1)*(h2m0./meffn.*kgrid.^2/qel);
    %
    lmbV=SBV(:,1);
    SBVtmp=lmbV*ones(1,nk)+ones(nvh,1)*(h2m0./meffp.*kgrid.^2/qel);
    
    Ban.SBC=SBCtmp;
    Ban.SBV=SBVtmp;
    
    % From 1990 Bava, Debernardi, Lin report, p. 21, eq. (39)
    % the cos(2*phi) terms can be ignored
    M2d=ti*R2.*(SHt); % for gain and spontaneous emission: no light holes
    M2esd=ti*R2.*(2*SHt)/3; % for spontaneous recombination: no light holes
    
    Ban.M2d=M2d(:,1)*ones(1,nk);
    Ban.M2esd=M2esd(:,1)*ones(1,nk);

    for indnk=1:nk
        Ban.SH(:,indk,:)=diag(squeeze(SH(:,1,:)));
        Ban.SL(:,indk,:)=diag(squeeze(SL(:,1,:)));
    end
    
else
    
    % From 1990 Bava, Debernardi, Lin report, p. 21, eq. (39)
    % the cos(2*phi) terms can be ignored
    M2d=ti*R2.*(SHt + SLt/3); % for gain and spontaneous emission
    M2esd=ti*R2.*(2*SHt + SLt)/3; % for spontaneous recombination

    Ban.SBC=SBC;
    Ban.SBV=SBV;
    
    Ban.M2d=M2d;
    Ban.M2esd=M2esd;
    
    Ban.SH=SH;
    Ban.SL=SL;
end

save(mesh.fileName,'Ban','mesh')



%==========================================================================
% Plot dispersion curves
%==========================================================================
if mode.iplot==1
    figure(2)
    set(gcf,'Position',[1060 512 560 420])
    hold on
    grid on
    plot(kgrid*1e-9,(-Ban.SBV),'LineWidth',2)
    set(gca,'FontSize',14,'FontName','Arial','Box','on')
    xlabel('k_{||}, nm^{-1}')
    ylabel('Energy, eV')
    title('Valence subbands')
    
    figure(3)
    set(gcf,'Position',[633 59 560 420])
    hold on
    grid on
    plot(kgrid*1e-9,(Ban.SBC),'LineWidth',2)
    set(gca,'FontSize',14,'FontName','Arial','Box','on')
    xlabel('k_{||}, nm^{-1}')
    ylabel('Energy, eV')
    title('Valence subbands')
end

