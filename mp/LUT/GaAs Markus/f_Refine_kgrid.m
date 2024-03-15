
function Ban=f_Refine_kgrid(Ban,mesh,mode)

% This function refines the kgrid, and then the conduction/valence
% subbands, in order to achieve better integration in k. Subbands are quite
% smooth, therefore spline interpolation is very effective for such
% purpose.

s_LoadConstants

kgrid=Ban.kgrid;
nc=1:mesh.ncb;
nv=1:mesh.nvb;
V0=mesh.V0;

SBC=Ban.SBC;
SBV=Ban.SBV;
M2d=Ban.M2d;
M2esd=Ban.M2esd;
num_kvectors=length(Ban.kgrid);
SH=Ban.SH;
nvh=Ban.nvh;

Mul_k=10; % Number of k points
if mode.ifit==0
    Mul_k=1; % Number of k points
end

kx1=linspace(kgrid(1),kgrid(end/2),num_kvectors/2*Mul_k); % k grid points along x, angstrom

kgridf=kx1; % k vectors grid, conversion to 1/m

kdkf=[0 kgridf(2:end).*diff(kgridf)]';
kdkf(end)=kdkf(end)/2;

Ban.kdkf=kdkf;
Ban.kgridf=kgridf;

for ke=nc
    SBCf(ke,:)=spline(kgrid,SBC(ke,:),kgridf);
end
Egc=SBCf.';

for kv=nv
    SBVf(kv,:)=spline(kgrid,SBV(kv,:),kgridf);
    for ke=nc
%         Egc(:,ke)=spline(kgrid,Ec(:,ke),kgridf);
        indb=kv+(ke-1)*length(nv);
        M2df(indb,:)=spline(kgrid,M2d(indb,:),kgridf);
        M2esdf(indb,:)=spline(kgrid,M2esd(indb,:),kgridf);
    end
end

Ban.SBCf=SBCf;
Ban.SBVf=SBVf;

for kv=nv
    indBound_k=find(SBVf(kv,:)>-V0);
    for ke=nc
        if length(indBound_k)>0
            indb=kv+(ke-1)*length(nv);
            M2df(indb,indBound_k)=0;
            M2esdf(indb,indBound_k)=0;
        end
    end
end

ECVf=[];
for ke=nc
    ECVf=[ECVf; SBVf+ones(nvh,1)*SBCf(ke,:)];
end

Ban.ECVf=ECVf;
Ban.M2df=M2df;
Ban.M2esdf=M2esdf;
