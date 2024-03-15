
function [lmbV,xvV]=solve_heavyVband_AlGaAs(mesh)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Finite element method for the solution of the 1-D stationary Schrhoedinger'sdffopksdopfsdf equation with
% Dirichlet boundary conditions.
%
% Alberto Tibaldi, 18/08/2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Physical constants
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
qel = 1.6021766208e-19;    % electron charge, C
c_light = 2.99792458e8;     % speed of light, m/s
kB = 1.3806488e-23;         % Boltzmann constant, J/K
h=6.626070040e-34;          % Planck constant, J*s
hbar = h/(2*pi);            % reduced Planck constant, J*s
mu0 = 4*pi*1e-7;            % magnetic permeability constant, H/m
eps0 = 1/(mu0*c_light^2);   % dielectric permittivity constant, F/m
m0 = 9.10938188e-31;        % Electron mass, kg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mesh indices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
in=1:mesh.nn; % node indexes
in1=in(1:end-1); in2=in(2:end); % "left" nodes for each element
inr=[in1 in1 in2 in2]; % assembling "rows"
inc=[in1 in2 in1 in2]; % assembling "columns"
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mesh parameters generation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
in=1:mesh.nn; % node indexes
in1=in(1:end-1); in2=in(2:end); % "left" nodes for each element
mesh.ne=mesh.nn-1; % number of mesh elements (nodes - 1)
mesh.z=linspace(0,mesh.L,mesh.nn); % mesh definition: [0,pi] interval
mesh.zc=(mesh.z(in1)+mesh.z(in2))/2; % Center points, m
mesh.Le = diff(mesh.z); % lengths of elements
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assembling matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Discretized Hamiltonian operator AlGaAs
% meff0=0.067;           % Electron effective DoS mass for xmol=0 (GaAs)
% meff1=0.124;            % Electron effective DoS mass for xmol=1 (AlAs)
% invmeff=((1-mesh.xmol)/meff0+mesh.xmol/meff1)/m0;
%mesh.xmol is zero in well so only GaAs!

% Discretized Hamiltonian operator InGaAs
%barrier
gamma1InAs=16.5;% Kim 2009
gamma2InAs=6.5;% Kim 2009
gamma1AlAs= 3.45; % Luttinger parameter, 1995Chuang, p. 709
gamma2AlAs= 0.68; % Luttinger p
gamma1GaAs = 6.85; % Luttinger parameter, 1995Chuang, p. 709
gamma2GaAs= 2.10; % Luttinger parameter, 1995Chuang, p. 709

gamma1well=mesh.xmol_well*gamma1GaAs+(1-mesh.xmol_well)*gamma1InAs;
gamma2well=mesh.xmol_well*gamma2GaAs+(1-mesh.xmol_well)*gamma2InAs;
gamma1barr=mesh.xmol_barrier*gamma1GaAs+(1-mesh.xmol_barrier)*gamma1AlAs;
gamma2barr=mesh.xmol_barrier*gamma2GaAs+(1-mesh.xmol_barrier)*gamma2AlAs;
mhwell=1/(gamma1well-2*gamma2well);
mlwell=1/(gamma1well+2*gamma2well);
mhbarr=1/(gamma1barr-2*gamma2barr);
mlbarr=1/(gamma1barr+2*gamma2barr);

invmeff=ones(1,mesh.ne)*(mlbarr*m0)^-1;
%well
vz_offset=0;
for indQW=1:length(vz_offset)
    ii=((mesh.xc>=(mesh.L/2-mesh.Lz/2+vz_offset(indQW))) & (mesh.xc<=(mesh.L/2+mesh.Lz/2+vz_offset(indQW))));
    invmeff(ii)=(mlwell*m0).^-1;
end

%  invmeff=1./(invmeff);
%
h2m0=hbar^2./(2.*qel).*(invmeff);
M11=h2m0./mesh.Le;
M12=-h2m0./mesh.Le;
M21=-h2m0./mesh.Le;
M22=h2m0./mesh.Le;
MM=[M11 M12 M21 M22];
Hmat=sparse(inr,inc,MM,mesh.nn,mesh.nn); % assembling Hamiltonian matrix
%
% Discretized potential
V=-mesh.evb2;
M11=2/6*V.*mesh.Le; % this is M_{n,n}
M12=1/6*V.*mesh.Le; % this is M_{n,n+1}
M21=1/6*V.*mesh.Le; % this is M_{n+1,n}
M22=2/6*V.*mesh.Le; % this is M_{n+1,n+1}
MM=[M11 M12 M21 M22];
Vmat=sparse(inr,inc,MM,mesh.nn,mesh.nn); % assembling potential matrix
%
% Mass matrix
M11=2/6*mesh.Le; % this is M_{n,n}
M12=1/6*mesh.Le; % this is M_{n,n+1}
M21=1/6*mesh.Le; % this is M_{n+1,n}
M22=2/6*mesh.Le; % this is M_{n+1,n+1}
MM=[M11 M12 M21 M22];
Mmat=sparse(inr,inc,MM,mesh.nn,mesh.nn); % assembling mass matrix
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solution of the eigenproblem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% eigs: sparse eigenproblem for 10 eigenvalues
opts.disp=0;
targetEigVl=0;
[xvV,lmbV]=eigs(Hmat+Vmat,Mmat,mesh.nvinit,targetEigVl,opts);
[lmbV,indEigVl]=sort(diag(real(lmbV))); % sorting eigenvalues
xvV=xvV(:,indEigVl);
% xvC'
