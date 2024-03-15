
function [Ban,mesh]=f_ComputeQWSubbands(mesh,mode,Qw_sel)

%==========================================================================
% Definition of FEM grid and structure parameters
%==========================================================================

vz_offset=mesh.vz_offset;
% vz_offset=[-mesh.Lz,mesh.Lz]; % 2 QWs case
%vz_offset=[(mesh.barr1QW+mesh.well1QW+mesh.barr2QW+mesh.well2QW+mesh.barr2QW),mesh.barr2QW+mesh.well3QW+mesh.barr3QW,mesh.barr3QW]; % 3 QWs case

s_LoadConstants
if Qw_sel==0
s_GaAsConstants
else
 s_InGaAsConstants
end
conv=h*c_light*10^9/qel; 
mesh.ne=mesh.nn-1; % Number of elements
mesh.x=linspace(0,mesh.L,mesh.nn); % Node coordinates, m
mesh.le=mesh.x(2:end)-mesh.x(1:end-1); % Edge length, m
mesh.xc=(mesh.x(1:end-1)+mesh.x(2:end))/2; % Center points, m
xmol_well=mesh.xmol_well;
xmol_barrier=mesh.xmol_barrier;
if Qw_sel==0
Delc=mesh.Qc*xmol_barrier;
Delh=(1-mesh.Qc)*xmol_barrier;

C0=Delc*(mesh.DeltaEg); % Conduction band offset, eV
V0=-Delh*(mesh.DeltaEg); % Valence band offset, eV

mesh.C0=C0;
mesh.V0=V0;
mesh.ecb=C0*ones(1,mesh.ne); % Conduction band offset vector
mesh.evb=V0*ones(1,mesh.ne); % Valence band offset vector
end

if Qw_sel==1
mesh.ecb=mesh.C0*ones(1,mesh.ne); % Conduction band offset vector
mesh.evb=mesh.V0*ones(1,mesh.ne); % Valence band offset vector    
C0=mesh.C0;
V0=mesh.V0;
end
mesh.xmol=xmol_barrier*ones(1,mesh.ne);
mesh.meffn=mesh.meffn_b*ones(1,mesh.ne)

%well(1)=mesh.well1QW;
%well(2)=mesh.well2QW;
%well(3)=mesh.well3QW;

%vz_offset=mesh.vz_offset;
well=mesh.well;

for indQW=1:length(vz_offset)
    ii=((mesh.xc>=(mesh.L-well(indQW)-vz_offset(indQW)) & (mesh.xc<=(mesh.L-vz_offset(indQW)))));
    mesh.evb(ii)=0;
    mesh.ecb(ii)=0;
    mesh.xmol(ii)=xmol_well;
    mesh.meffn(ii)=mesh.meffn_w;
end

if mode.iplot==1
    %
    figure(1)
%     subplot(2,1,1)
    axis on
    grid on
    hold on
    box on
    plot(mesh.xc*1e9,(mesh.ecb+mesh.Eg),'b.')
    plot(mesh.xc*1e9,mesh.evb,'r.')
%       plot(mesh.xc*1e9,mesh.ecb+mesh.Eg,'b','LineWidth',2)
%     plot(mesh.xc*1e9,mesh.evb,'r','LineWidth',2)
    xlabel('z (nm)')
    ylabel('Valence band structure, eV')
    xlim([0 mesh.L*1e9])
%     subplot(2,1,2)
%         axis onb
%         grid on
%         hold on
%         box on
%         plot(mesh.xc*1e9,mesh.meffn,'bo')
%         xlabel('z (nm)')
%        ylabel('effective mass')
%       xlim([0 mesh.L*1e9])

end
% a=mesh.ecb;
% b=mesh.evb;
% c=mesh.Eg;
% save cb a 
% save vb b
% save gape c
%  load cb
% load vb
% load gape
% mesh.ecb=a;
% mesh.evb=b;
% mesh.Eg=c;
% plot(mesh.xc*1e9,mesh.ecb+mesh.Eg,'g--','LineWidth',2)
%     plot(mesh.xc*1e9,mesh.evb,'k--','LineWidth',2)
% legend('Cb strain','Vb strain','Cb no str','Vb no str')
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
if Qw_sel==0
mesh.targetC=min(mesh.ecb)+1.5; % target for eigenvalue solver, conduction band
[lmbVtmp,xvVtmp,lmbCtmp,xvCtmp]=solve_kp44_AlGaAs(mesh,0,0,mode);
[lmbCtmps,xvCtmps]=solve_Cband_AlGaAs(mesh);
end
if Qw_sel==1
  mesh.targetC=min(mesh.ecb)+mesh.Eg+0.01;  
  [lmbVtmp,xvVtmp,lmbCtmp,xvCtmp,mesh]=solve_kp44_InGaAs(mesh,0,0,mode);
  mesh.evb1=mesh.evb;
  mesh.evb2=mesh.evb;
  mesh.ecb1=mesh.ecb+mesh.Eg;
  for indQW=1:length(vz_offset)
      ii=((mesh.xc>=(mesh.L-well(indQW)-vz_offset(indQW)) & (mesh.xc<=(mesh.L-vz_offset(indQW)))));
      mesh.evb1(ii)=mesh.v1;
      mesh.evb2(ii)=mesh.v2;
      mesh.ecb1(ii)=mesh.c1;
      mesh.xmol(ii)=xmol_well;
      mesh.meffn(ii)=mesh.meffn_w;
  end
%   
%  for indQW=1:length(vz_offset)
%        ii=~((mesh.xc>=(mesh.L-well(indQW)-vz_offset(indQW)) & (mesh.xc<=(mesh.L-vz_offset(indQW)))));
%       mesh.evb1(ii)=mesh.v1p+mesh.evb(ii);
%       mesh.evb2(ii)=mesh.v2p+mesh.evb(ii);
%    %mesh.ecb1(ii)=mesh.c1p; % not to be included, not considtent with the
%     %model of Delta Ec and Delta Ev
% end
  [lmbCtmps,xvCtmps]=solve_Cband_AlGaAs(mesh);
end

%% discard bands not confined, change of target ncivit and nvit to avoid their computation
 maxCb=mesh.ecb1(1);
 minVb=mesh.evb2(1);
 p1cb=find(lmbCtmp<maxCb);
 p2vb=find(lmbVtmp>minVb+0.02); 
 
mc=p1cb(end); % number of conduction subbands to be computed
mv=p2vb(end); % number of valence subbands to be computed
%recompute starting point

xvC=reshape(xvCtmp,4,mesh.nn,mesh.ncinit);
xvV=reshape(xvVtmp,4,mesh.nn,mesh.nvinit);
txt={'first CB','second CB','third CB','fourth CB'};
% 
% hfig=figure
% set(hfig,'pos',[       32          10        1867         974])
for indLivello=1:mc
subplot(1,mc,indLivello)
xv = reshape(xvCtmp(:,indLivello),4,mesh.nn);
plot(mesh.x*1e9,abs(sum(xv,1).^2),'linewidth',1.2)
hold on
plot(mesh.x*1e9,abs(xvCtmps(:,indLivello).^2),'k--','linewidth',1.2)
legend(['Kp, \lambda =' num2str(conv/lmbCtmp(indLivello),4)],['Schrodinger, \lambda =' num2str(conv/lmbCtmps(indLivello),4)]);
%title(txt{indLivello});
xlabel('z (nm)')
ylabel('{|\Psi|}^2')
box on
grid on
set(gca,'FontSize',12,'FontName','Arial','box','on')
xlim([0 mesh.L*1e9])
end
 
 
 
 
 
 
figure
 
 for indLivello=1:mc
 %figure(hfig1)
 xv = reshape(xvCtmp(:,indLivello),4,mesh.nn);
% P=sum(xv,1);
 P=xv(1,:); 
 plot(mesh.x*1e9,P*sign(P(20)),'linewidth',1.2)
hold on
 P=xvCtmps(:,indLivello);
 plot(mesh.x*1e9,P*sign(P(20)),'k--','linewidth',1.2)
 legend(['Kp, \lambda =' num2str(conv/lmbCtmp(indLivello),4)],['Schrodinger, \lambda =' num2str(conv/lmbCtmps(indLivello),4)]);
 %title(txt{indLivello});
 xlabel('z (nm)')
 ylabel('\Psi')
 box on
 grid on
 set(gca,'FontSize',12,'FontName','Arial','box','on')
 xlim([0 mesh.L*1e9])
 
 
end
 
 
 pausak
 
figure
hold on
 axis on
    grid on
    hold on
    box on
    plot(mesh.xc*1e9,(mesh.ecb1),'b.')
    plot(mesh.xc*1e9,mesh.evb1,'r.')
    plot(mesh.xc*1e9,mesh.evb2,'g.')

for indLivello=1:length(lmbCtmp)
plot(mesh.x*1e9,ones(size(mesh.x))*lmbCtmp(indLivello),'--')
end
for indLivello=1:length(lmbVtmp)
plot(mesh.x*1e9,ones(size(mesh.x))*lmbVtmp(indLivello),'--')
end
    xlabel('z (nm)')
    ylabel('band structure, eV')
    xlim([0 mesh.L*1e9])
xlim([0 mesh.L*1e9])

pausak

xvV=reshape(xvVtmp,4,mesh.nn,mesh.nvinit);

%numliv=3; no strain
if Qw_sel==1
  [lmbVHtmps,xvVHtmps]=solve_heavyVband_AlGaAs(mesh);
  [lmbVLtmps,xvVLtmps]=solve_ligVband_AlGaAs(mesh);

lmbVtmpst=[lmbVHtmps;lmbVLtmps];
xvVtmpst=[xvVHtmps xvVLtmps]; 
[ lmbVtmps, ord]=sort(lmbVtmpst);
xvVtmps=xvVtmpst(:,ord);

hfig=figure
set(hfig,'pos',[       32          10        1867         974])

for indLivello=1:mv
subplot(round(mv/2),round(mv/2),indLivello);

xv = reshape(xvVtmp(:,indLivello),4,mesh.nn);
plot(mesh.x*1e9,abs(sum(xv,1).^2),'linewidth',1.2)
hold on
plot(mesh.x*1e9,xv(3,:).^2,'--','linewidth',1.2)

plot(mesh.x*1e9,abs(xvVtmps(:,indLivello).^2),'k--','linewidth',1.2)
legend(['Kp, \lambda =' num2str(lmbVtmp(indLivello),4)],['Schrodinger, \lambda =' num2str(-lmbVtmps(indLivello),4)]);
legend('Location','best')

hold on
xlabel('z (nm)')
ylabel('{|\Psi|}^2')
box on
grid on
set(gca,'FontSize',12,'FontName','Arial','box','on')
xlim([0 mesh.L*1e9])
end

pausak




hfig=figure
set(hfig,'pos',[       32          10        1867         974])

for indLivello=1:mv
subplot(round(mv/2),round(mv/2),indLivello);

xv = reshape(xvVtmp(:,indLivello),4,mesh.nn);
P=sum(xv,1);
plot(mesh.x*1e9,abs(P).^2,'linewidth',1.2)
hold on
P=xvVtmps(:,indLivello);
plot(mesh.x*1e9,abs(P).^2,'k--','linewidth',1.2)
legend(['Kp, \lambda =' num2str(lmbVtmp(indLivello),4)],['Schrodinger, \lambda =' num2str(-lmbVtmps(indLivello),4)]);
legend('Location','best')

hold on
xlabel('z (nm)')
ylabel('{|\Psi|}^2')
box on
grid on
set(gca,'FontSize',12,'FontName','Arial','box','on')
xlim([0 mesh.L*1e9])
end

%keyboard
hfig=figure
set(hfig,'pos',[       32          10        1867         974])

for indLivello=1:mv
subplot(round(mv/2),round(mv/2),indLivello);

xv = reshape(xvVtmp(:,indLivello),4,mesh.nn);
P=xv(2,:);
%P=sum(xv,1);
plot(mesh.x*1e9,P*sign(P(1)),'linewidth',1.2)
hold on
P1=xv(3,:);
%P=sum(xv,1);
plot(mesh.x*1e9,P1*sign(P1(1)),'linewidth',1.2)
P=xvVtmps(:,indLivello);
plot(mesh.x*1e9,P*sign(P(1)),'k--','linewidth',1.2)
legend(['Kp, \lambda =' num2str(lmbVtmp(indLivello),4)],['Schrodinger, \lambda =' num2str(-lmbVtmps(indLivello),4)]);
legend('Location','best')

hold on
xlabel('z (nm)')
ylabel('{\Psi}')
box on
grid on
set(gca,'FontSize',12,'FontName','Arial','box','on')
xlim([0 mesh.L*1e9])
end


end
% con strain si vede che le psi non sono simmetriche quando ha a che fare l-accopiamento con elettroni, quindi banda leggera e tutte cb, mentre vb pesante no infatti e simmetrica
%%keyboard
pausak


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

keyboard




disp('Computing subbands')

for ik=1:mesh.num_kvectors
    
    fprintf('k, 1/A %e %e\n',[kx(ik) ky(ik)])
    
    if Qw_sel==0
    [lmbV,xvV,lmbC,xvC]=solve_kp44_AlGaAs(mesh,kx(ik),ky(ik),mode);
    end
    
    if Qw_sel==1
  [lmbV,xvV,lmbC,xvC]=solve_kp44_InGaAs(mesh,kx(ik),ky(ik),mode);
end
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

%indBound_V=find(SBV(:,1)<=-minVb-0.02);
%indBound_C=find(SBC(:,1)<=maxCb);
%mesh.ncb=length(indBound_C);
%mesh.nvb=length(indBound_V);

indBound_V=find(SBV(:,1)<=-mesh.evb2(1)*.9);
indBound_C=find(SBC(:,1)<=mesh.ecb1(1)-(max(mesh.ecb1)-min(mesh.ecb1))*.1);
mesh.ncb=length(indBound_C);
mesh.nvb=length(indBound_V);



%XVV=XVV(:,indBound_V,:);
%SBV=SBV(indBound_V,:);
%XVC=XVC(:,indBound_C,:);
%SBC=SBC(indBound_C,:);

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

XV1old=XV1;
XV2old=XV2;
SBVold=SBV;

h=figure, 
set(h,'pos',[346         453        1308         525])
subplot(121)
plot(kgrid*1e-9,SBC)
subplot(122)
plot(kgrid*1e-9,SBV)


'vedi bande prima di riduzione', keyboard

XV1=XV1(:,indBound_V,:);
XV2=XV2(:,indBound_V,:);
%XVV=XVV(:,indBound_V,:);
SBV=SBV(indBound_V,:);
XVC=XVC(:,indBound_C,:);
SBC=SBC(indBound_C,:);

h=figure, 
set(h,'pos',[346         453        1308         525])
subplot(121)
plot(kgrid*1e-9,SBC)
subplot(122)
plot(kgrid*1e-9,SBV)
pausak

%======================================================================
% Computing dipole momentum
%======================================================================
Eg=qel*mesh.Eg; % bandgap, J
Delta=qel*mesh.Delta; % spin-orbit coupling, J
mc=mesh.meffn_w*m0; % conduction band effective mass, kg

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

XV1s=XV1;
XV2s=XV2;
XVCs=XVC;

TipoH=max(abs(squeeze(XV1s(:,:,1))))-max(abs(squeeze(XV2s(:,:,1))));
iHH=find(TipoH>0);
iLH=find(TipoH<0);



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

if nc==1
 SHp=SH;
 SLp=SL;
else
 SHp=[];
 SLp=[];
 for se=1:nc
  SHp=[SHp; SH(:,:,se)];
  SLp=[SLp; SL(:,:,se)];
 end 
end

figure, plot(kgrid*1e-9,abs(SHp.^2)), hold on
ax = gca;
ax.ColorOrderIndex = 1;
plot(kgrid*1e-9,abs(SLp.^2),'--'),
xlabel('k_t 1/nm')
ylabel('Overlap HH and LH (--)')

icontr=0;
if icontr==1
R1=reshape(XV1(:,1:3:end),661,22);
R2=reshape(XV1(:,2:3:end),661,22);
R3=reshape(XV1(:,3:3:end),661,22);
L3=reshape(XV2(:,3:3:end),661,22);
L2=reshape(XV2(:,2:3:end),661,22);
L1=reshape(XV2(:,1:3:end),661,22);
figure, plot(mesh.x*1e9,R1(:,1:7:end))
hold on
ax = gca;
ax.ColorOrderIndex = 1; plot(mesh.x*1e9,L1(:,1:7:end),'.')
xlim([25 40])
pausak
figure, plot(mesh.x*1e9,R2(:,1:7:end))
hold on
ax = gca;
ax.ColorOrderIndex = 1; plot(mesh.x*1e9,L2(:,1:7:end),'.')
xlim([25 40])
title('L LH1: --- HHpart, ... LH part')
pausak
figure, plot(mesh.x*1e9,R3(:,1:7:end))
hold on
ax = gca;
ax.ColorOrderIndex = 1; plot(mesh.x*1e9,L3(:,1:7:end),'.')
xlim([25 40])
title('U HH2: --- HHpart, ... LH part')
end

'fine overlap', keyboard

sih=size(SH);
SHt=[]; SLt=[]; ECV=[];

for ke=1:nc
    SHt=[SHt; SH(:,:,ke).^2];
    SLt=[SLt; SL(:,:,ke).^2];
    ECV=[ECV; SBV+ones(sih(1),1)*SBC(ke,:)];
end

% 
% SHd=SHt;
% SLd=SLt;
% save val SHd SLd;
%======================================================================
% Computing potentials for many-body effects. Pote e linhard solo per k=0
%======================================================================
XVC=reshape(XVC,mesh.nn,mesh.ncb,nk);
XVC=squeeze(XVC(:,:,1));
XV1d=reshape(XV1,mesh.nn,mesh.nvb,nk); % per plot a k diversi
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

for ice = 1:mesh.ncb
    fe=XVC(:,ice);
    PPEcV{ice} = spline(kgrid,SBC(ice,:));

for ic = 1:mesh.ncb
    fi=XVC(:,ic);
    f14 = (fe.*fi).';
    for iq = 1:length(kgrid)
        qD = kgrid(iq) + qs;
        % Unscreened conduction potential, Veprek PhD, (5.2)
        Vcc(iq) = qel^2/(2*eps0*ER_STA*qD)*(dz^2)*(f14*exp(-qD*XY)*f14.');
        if iq==length(kgrid)*10
         'iq ', iq, ice, ic, 
         PIN=(diag(f14)*exp(-qD*XY)*diag(f14));
         figure, surf(PIN), shading flat, view(2), pausak
        end 
    end
    % Spline coefficients for conduction potential
%    PPVccV{ic,ice} = spline(kgrid,Vcc);
    VS(:,ic,ice)=Vcc;
end
 Vcs=sum(VS(:,:,ice),2);
 VSE(:,ice)=Vcs;
'ice', ice, pausak
    PPVccV{ice} = spline(kgrid,Vcs);
end

for ive = iHH
%for ive = 4:6
    PPEvV{ive} = spline(kgrid,SBV(ive,:));
    fe=XVV(:,ive);
for iv = iHH
    fi=XVV(:,iv);
    f14 = (fe.*fi).';
    for iq = 1:length(kgrid)
        qD = kgrid(iq) + qs;
        % Unscreened conduction potential, Veprek PhD, (5.2)
        Vcc(iq) = qel^2/(2*eps0*ER_STA*qD)*(dz^2)*(f14*exp(-qD*XY)*f14.');
        if iq==length(kgrid)*10
         'iq ', iq, ive, iv, 
         PIN=(diag(f14)*exp(-qD*XY)*diag(f14));
         figure, surf(PIN), shading flat, view(2), pausak
        end 
    end
    % Spline coefficients for conduction potential
%    PPVccV{ic,ice} = spline(kgrid,Vcc);
    VSh(:,iv,ive)=Vcc;
end
 Vcs=sum(VSh(:,:,ive),2);
 VSH(:,ive)=Vcs;
%'ive', ive, pausak
    PPVvvV{ive} = spline(kgrid,Vcs);
end


for ive = iLH
    PPEvV{ive} = spline(kgrid,SBV(ive,:));
    fe=XVV(:,ive);
for iv = iLH
    fi=XVV(:,iv);
    f14 = (fe.*fi).';
    for iq = 1:length(kgrid)
        qD = kgrid(iq) + qs;
        % Unscreened conduction potential, Veprek PhD, (5.2)
        Vcc(iq) = qel^2/(2*eps0*ER_STA*qD)*(dz^2)*(f14*exp(-qD*XY)*f14.');
        if iq==length(kgrid)*10
         'iq ', iq, ice, ic, 
         PIN=(diag(f14)*exp(-qD*XY)*diag(f14));
         figure, surf(PIN), shading flat, view(2), pausak
    end
    end
    % Spline coefficients for conduction potential
%    PPVccV{ic,ice} = spline(kgrid,Vcc);
    VSh(:,iv,ive)=Vcc;
end
 Vcs=sum(VSh(:,:,ive),2);
 VSH(:,ive)=Vcs;
%'ice', ice, pausak
    PPVvvV{ive} = spline(kgrid,Vcs);
end

figure, semilogy(kgrid*1e-9,VSE)
hold on
semilogy(kgrid*1e-9,VSH,'--')
pausak

figure, subplot(211)
plot((mesh.x-mesh.x(ceil(end/2)))*1e9,XVC), subplot(212)
plot((mesh.x-mesh.x(ceil(end/2)))*1e9,sum(XVC.^2,2)), pausak

figure, subplot(211)
plot((mesh.x-mesh.x(ceil(end/2)))*1e9,XVV), subplot(212)
plot((mesh.x-mesh.x(ceil(end/2)))*1e9,sum(XVV.^2,2)), pausak

%for iv = 1:mesh.nvb
%    f23 = abs(XVV(:,iv).').^2;
%    PPEvV{iv} = spline(kgrid,SBV(iv,:));
%    for iq = 1:length(kgrid)
%        qD = kgrid(iq) + qs;
%        % Unscreened valence potential, Veprek PhD, (5.2)
%        Vvv(iq) = qel^2/(2*eps0*ER_STA*qD)*(dz^2)*(f23*exp(-qD*XY)*f23.');
%    end
%    % Spline coefficients for valence potential
%    PPVvvV{iv} = spline(kgrid,Vvv);
%end

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


'qui dipolo 3QW;', keyboard
save(mesh.fileName,'Ban','mesh')
% 
% load strainCB 
% load strainVB

%==========================================================================
% Plot dispersion curves
%==========================================================================
if mode.iplot==1
    figure
 
    hold on
    grid on
    p1=plot(kgrid*1e-9,(-Ban.SBV),'LineWidth',2)
    %p1=plot(kgrid*1e-9,(-Ban.SBV),'k','LineWidth',2)
    hold on
   % p2=plot(kgrid*1e-9,(-a),'r--','LineWidth',2)
    set(gca,'FontSize',14,'FontName','Arial','Box','on')
    xlabel('k_{||}, nm^{-1}')
    ylabel('Energy, eV')
    title('Valence subbands')
    %legend([p1(1) p2(1)],'no strain','strain')
      figure
 
    hold on
    grid on
    %p3=plot(kgrid*1e-9,(Ban.SBC),'k','LineWidth',2)
    p3=plot(kgrid*1e-9,(Ban.SBC),'LineWidth',2)
    hold on
   % p4=plot(kgrid*1e-9,(b),'r--','LineWidth',2)
    set(gca,'FontSize',14,'FontName','Arial','Box','on')
    xlabel('k_{||}, nm^{-1}')
    ylabel('Energy, eV')
    title('Conduction subbands')
    %legend([p3(1) p4(1)],'no strain','strain')
end
% a=Ban.SBV;
% save strainVB a
%  b=Ban.SBC;
% save strainCB b
