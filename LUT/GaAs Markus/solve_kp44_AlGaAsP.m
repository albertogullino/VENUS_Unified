function [lmbV,xvV,lmbC,xvC]=solve_kp44_AlGaAs(mesh,kx,ky)
%
p_asy=1; % Burt-Foreman ordering
% p_asy=0.5; % Symmetric ordering
%
nB = 4;
nn=mesh.nn; % number of mesh nodes
%
I = eye(nB);
%
Amat=spalloc(nn*nB,nn*nB,nB^2*(3*nn-2));
Bmat=spalloc(nn*nB,nn*nB,nB*(3*nn-2));
%
for ie=1:mesh.ne % Loop over mesh elements
    le=mesh.le(ie); % element length, m
    % bulk Hamiltonian
    [H0,H1l,H1r,H2]=assem_kp44_AlGaAs(kx,ky,mesh.xmol(ie));
    %
    H0 = H0 + mesh.evb(ie)*I;
    H1L=H1l*p_asy+H1r*(1-p_asy);
    H1R=H1r*p_asy+H1l*(1-p_asy);
    %
    NN11=1/3*le;
    NN12=1/6*le;
    NN22=1/3*le;
    NxNx11=1/le;
    NxNx12=-1/le;
    NxNx22=1/le;
    NxN=1i/2*[-1 -1; 1  1];
    NNx=NxN';
    %
    i1=(1:nB)+(ie-1)*nB;
    i2=(1:nB)+ie*nB;
    %
    % assem
    Amat(i1,i1)=Amat(i1,i1)+H2*NxNx11;
    Amat(i1,i2)=Amat(i1,i2)+H2*NxNx12;
    Amat(i2,i1)=Amat(i2,i1)+H2*NxNx12;
    Amat(i2,i2)=Amat(i2,i2)+H2*NxNx22;
    %
    Amat(i1,i1)=Amat(i1,i1)+H1L*NxN(1,1);
    Amat(i1,i2)=Amat(i1,i2)+H1L*NxN(1,2);
    Amat(i2,i1)=Amat(i2,i1)+H1L*NxN(2,1);
    Amat(i2,i2)=Amat(i2,i2)+H1L*NxN(2,2);
    %
    Amat(i1,i1)=Amat(i1,i1)+H1R*NNx(1,1);
    Amat(i1,i2)=Amat(i1,i2)+H1R*NNx(1,2);
    Amat(i2,i1)=Amat(i2,i1)+H1R*NNx(2,1);
    Amat(i2,i2)=Amat(i2,i2)+H1R*NNx(2,2);
    %
    Amat(i1,i1)=Amat(i1,i1)+(H0)*NN11;
    Amat(i1,i2)=Amat(i1,i2)+(H0)*NN12;
    Amat(i2,i1)=Amat(i2,i1)+(H0)*NN12;
    Amat(i2,i2)=Amat(i2,i2)+(H0)*NN22;
    %
    Bmat(i1,i1)=Bmat(i1,i1)+I*NN11;
    Bmat(i1,i2)=Bmat(i1,i2)+I*NN12;
    Bmat(i2,i1)=Bmat(i2,i1)+I*NN12;
    Bmat(i2,i2)=Bmat(i2,i2)+I*NN22;
end

opts.disp=0; % verbose mode of the iterative eigensolver
targetC=mesh.targetC;
targetV=mesh.targetV;
%
[xvC,lmbC,flag]=eigs(Amat,Bmat,mesh.ncinit,targetC,opts);
if(flag~=0),disp('Warning! Eigensolver not converged!'),end
[lmbC,ind]=sort(diag(real(lmbC)),'ascend');
xvC=xvC(:,ind);
%
[xvV,lmbV,flag]=eigs(Amat,Bmat,mesh.nvinit,targetV,opts);
if(flag~=0),disp('Warning! Eigensolver not converged!'),end
[lmbV,ind]=sort(diag(real(lmbV)),'descend');
xvV=xvV(:,ind);
