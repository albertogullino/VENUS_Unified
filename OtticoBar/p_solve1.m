function mode=p_solve1(geom,mesh,mode)
%
% =============================================================================================100
% History: created Francesco Bertazzi
%          working Alberto Tibaldi
%
% Last update: 02 March 2017, still in writing
%
% Description: Poisson-Boltzmann/Fermi solver
%
% Graded heterointerfaces added, May 12, 2008
% Treatment of abrupt heterointerfaces, Nov 25, 2008
% Use of Fermi-Dirac statistics, Sept. 2016
% Incomplete ionization mode, Opt. 2016
% Quasi-3D (cylindrical) simulation, Nov. 2016
% Temperature-dependent material parameters, 01 Feb. 2017
%
% input variables:
% geom                     structure, geometry data created by FEM
% mesh                     structure, mesh data, see function RECTMESH
% mode.v0_ls               array(nm,nv): applied external voltages  
% mode.maxiter             maximum number of Newton iterations
% mode.maxiter_t           maximum number of Newton inner iterations
% mode.tolconv             required tolerance
% mode.report              flag: if true select verbose mode 
% mode.Zmat                impedance matrix of the linear network 
%                          embedding the nonlinear device
%
% output variables:
% mode.i_ls                array(1*nm): LS total current at last iteration, A/cm
% mode.v_ls                array(1*nm): LS voltage at last iteration, V
% mode.vv_ls               array(nm,nv): LS voltage, V 
% mode.ii_ls               array(nm,nv): LS total current, A/cm
% mode.phi                 array(1**nn): LS electrostatic potential, V
% mode.elec                array(1**nn): LS electron concentration, 1/cm^3
% mode.hole                array(1*nn): LS hole concentration, 1/cm^3
% mode.efieldx             array(1*nt): electric field alog x, V/cm
% mode.efieldy             array(1*nt): electric field alog y, V/cm
% mode.iconv               array(nv): 1 convergence achieved, 2 convergence failure 
% mode.res                 array(nv): residual error
% mode.Vt                  thermal voltage (potential scaling factor), V
%
% =============================================================================================100
% The residual vector is computed as: Kmat0*uvet+tvet=rvet
% The update vector is computed as: Jmat0*uvet=-rvet
% =============================================================================================100
s_LoadConstants
Vt=kB.*mesh.T./qel; % electrical equivalent in temperature, V
%
elapsed_time=cputime;
% =============================================================================================100 
% On first call, the routine LS_SOLVE finds the Poisson-Boltzmann (-Fermi) solution.
% LS_SOLVE may be restarted from a previous LS or DC solution; 
% Notice that the DC solution may be obtained by setting the number of harmonics to zero.
% neutrality solution =========================================================================100
if(not(isfield(mode,'phi'))) % this part of the code is executed only on first call 
%
mode.elec_neutr=zeros(1,mesh.nn); 
mode.phi_neutr =zeros(1,mesh.nn); 
mode.hole_neutr=zeros(1,mesh.nn);
ii = mesh.dop>0; jj = mesh.dop<0;
%
mode.elec_neutr(ii) = 1./2*(mesh.dop(ii)+sqrt(mesh.dop(ii).^2+4*mesh.nint(ii).^2)); 
mode.hole_neutr(ii) = mesh.nint(ii).^2./mode.elec_neutr(ii);
mode.hole_neutr(jj) = 1./2*(-mesh.dop(jj)+sqrt(mesh.dop(jj).^2+4*mesh.nint(jj).^2)); 
mode.elec_neutr(jj) = mesh.nint(jj).^2./mode.hole_neutr(jj);
mode.phi_neutr(ii)  = Vt(ii).*log(mode.elec_neutr(ii)./mesh.nint(ii)) + mesh.phi_r(ii);
mode.phi_neutr(jj)  = Vt(jj).*log(mode.elec_neutr(jj)./mesh.nint(jj)) + mesh.phi_r(jj);
mode.phi  = mode.phi_neutr; 
mode.elec = mode.elec_neutr; 
mode.hole = mode.hole_neutr;
%
mode.Vt=Vt;  
Vt_tr=pdeintrp(mesh.node,mesh.triangle(1:4,:),Vt.'); % T on triangles
mode.Vt_tr=Vt_tr;  
%
% Correct neutrality conditions for Fermi-Dirac and/or Incomplete
% ionization
[mode] = f_ComputeNeutrality(geom,mesh,mode);
end
% =============================================================================================100
% solution vector uvet: {phi,elec,hole,elec_t,v_ls,i_ls}
% vector indexes:
%
%    in    -> pointer to Poisson equation
%    in+qq -> pointer to voltage equations
%    in+rr -> pointer to current equations
%
nm=geom.nm;               % number of active contacts
%
nn=mesh.nn;               % -> elec eqs.
qq=nn;                    % -> circuit eqs.
rr=qq+nm;                 % -> current eqs.
%
% total number of equations/unknowns
neq=nn+2*nm;
%
% =============================================================================================100
% We first compute the solution vector with all the information we have at this point; 
% notice that here the solution is in frequency domain
uvet=zeros(1,neq); 
uvet(1:nn)=mode.phi;
% =============================================================================================100
%
%
% initializations
mode.vv_dd=[]; mode.ii_dd=[]; mode.Pst=0; 
mode.iconv=0;
%    
fprintf('\n Poisson solver\n')
v0_ls = 0;
mode.iconv = 0; % flag, 1 convergence, 2 convergence failure
iter = 0; % number of outer iterations 
iter_t = 0; % number of inner iterations
mode.res = Inf; % residual error
t = 1; % damping parameter
%
while(iter<mode.maxiter && iter_t<mode.maxiter_t), % outer while-loop
% assem jacobian matrix in time domain
[Kmat0,Jmat0,Jmat2,uvet,rvet,mode]=assem(geom,mesh,mode,uvet,v0_ls);
% solve matrix equation 
Jmat = Kmat0 + Jmat0 + Jmat2;
[R,C] = dgsequ(Jmat);
res = norm(R*rvet);
%==============================================================================================100
if(isnan(res)), error('residual error is NaN'), end
if(iter==0 && mode.report), fprintf('Iteration     Residual\n')
                            fprintf(sprintf('%4i%20.4e\n',iter,res)), end
if(iter>=1 && mode.report), fprintf('Iteration     Residual            t\n')
                            fprintf(sprintf('%4i%20.4e%15.2e\n',iter,res,t)), end
mode.res=res; % save residual
if(res<mode.tolconv), mode.iconv=1; mode.res = res; break, end
iter_t=0; t=1; % restore damping parameter
%=========================================
% solve matrix equation
Jmat = R*Jmat*C;
[Lmat,Umat,pm,qm] = lu(Jmat);
deltau = - C*(qm*(Umat\(Lmat\(pm*(R*rvet)))));
iter = iter + 1;
uvet = uvet + t*deltau.'; % update
end % outer while-loop
if(not(mode.iconv))
    fprintf('Warning: convergence failure\n')
end
%
% =============================================================================================100
% save LS solution
% =============================================================================================100
mode.phi=uvet(1:nn); 
%
mode.v_dd=uvet(qq+1); % saving voltage
mode.i_dd=uvet(rr+1); % saving current
mode.rvet=rvet.'; % save residual, debug only
fprintf('Elapsed time: %g, s\n',(cputime-elapsed_time))
mode.elapsed_time=(cputime-elapsed_time)/60;
mode.firstrun=0;
%**********************************************************************************************100
%
%
%
%**********************************************************************************************100
function [Kmat0,Jmat0,Jmat2,uvet,rvet,mode,tvet]=assem(geom,mesh,mode,uvet,v0_ls)
%
% define constants
s_LoadConstants
Vt=mode.Vt; % electrical equivalent in temperature, V
% =============================================================================================100
% Loading solution database
% =============================================================================================100
nt=mesh.nt;               % number of triangles
nm=geom.nm;               % number of active contacts 
nc=geom.nc;               % number of contacts 
%
nn=mesh.nn;               % -> elec eqs.
qq=nn;                    % -> circuit eqs.
rr=qq+nm;                 % -> current eqs.
%
% total number of equations/unknowns
neq=nn+2*nm;
%
% =============================================================================================100
% Loading and computing geometry parameters
% =============================================================================================100
% geometrical data
triangle=mesh.triangle; node=mesh.node; contact=mesh.contact; iq=mesh.iq;
% indexes to the nodes 1, 2, 3 of each triangle
in1=triangle(1,:); in2=triangle(2,:); in3=triangle(3,:); 
%
% x and y positions of node triangles
x1=node(1,in1); x2=node(1,in2); x3=node(1,in3);
y1=node(2,in1); y2=node(2,in2); y3=node(2,in3);
% distances from node triangles (FEM computations)
r3x=x2-x1; r1x=x3-x2; r2x=x1-x3;
r3y=-y2+y1; r1y=-y3+y2; r2y=-y1+y3;
% FEM basis functions
f23xx=r2x.*r3x; 
f23yy=r2y.*r3y; 
f23xy=r2x.*r3y; 
f32xy=r3x.*r2y;
% triangle side lengths
l1=sqrt(r1x.^2+r1y.^2); l2=sqrt(r2x.^2+r2y.^2); l3=sqrt(r3x.^2+r3y.^2);
xm1=1/2*(x2+x3); xm2=1/2*(x3+x1); xm3=1/2*(x1+x2);
ym1=1/2*(y2+y3); ym2=1/2*(y3+y1); ym3=1/2*(y1+y2);
% split points
xcc=1./2.*(x2.*f32xy+x1.*f32xy-y2.*f23yy-x3.*f23xy-x1.*f23xy+y3.*f23yy)./(f32xy-f23xy);
ycc=1./2.*(x2.*f23xx-y2.*f23xy-y1.*f23xy-x3.*f23xx+y3.*f32xy+y1.*f32xy)./(f32xy-f23xy);
% triangle mid-points
s1=sqrt((xm1-xcc).^2+(ym1-ycc).^2)*0; 
s2=sqrt((xm2-xcc).^2+(ym2-ycc).^2); 
s3=sqrt((xm3-xcc).^2+(ym3-ycc).^2); 
Se=abs(r2x.*r1y-r2y.*r1x)./2; % triangle area
% split areas (for each triangle, the area associated to node 1, 2, 3)
Se1=1/2.*(xm3.*ycc-xcc.*ym3+x1.*ym3-y1.*xm3+xcc.*ym2-xm2.*ycc-x1.*ym2+y1.*xm2);
Se2=1/2.*(xm1.*ycc-xcc.*ym1+x2.*ym1-y2.*xm1+xcc.*ym3-xm3.*ycc-x2.*ym3+y2.*xm3);
Se3=1/2.*(xm2.*ycc-xcc.*ym2+x3.*ym2-y3.*xm2+xcc.*ym1-xm1.*ycc-x3.*ym1+y3.*xm1);
%
% marking triangles outside semiconductor regions
semiconductor=find(geom.semiconductor);
it=not(ismember(triangle(4,:),semiconductor));
% by setting to zero Se1, Se2, Se3, we ensure that continuity equations 
% and trap equations are assembled only in the semiconductor regions! 
Se1(it)=0; Se2(it)=0; Se3(it)=0;
%
% =============================================================================================100
% Cylindrical coordinates correction
% =============================================================================================100
G1 = 1; G2 = 1; G3 = 1;
if((isfield(mode,'symmetry'))&&(strcmp(mode.symmetry,'Cylindrical-Y')))
xc = (x1+x2+x3)/3;
Se1=2*pi*xc.*Se1;
Se2=2*pi*xc.*Se2;
Se3=2*pi*xc.*Se3;
G1 = 2*pi*xc.*Gji(x2,x3,mode.tolGji); 
G2 = 2*pi*xc.*Gji(x3,x1,mode.tolGji); 
G3 = 2*pi*xc.*Gji(x1,x2,mode.tolGji);  end
%
if((isfield(mode,'symmetry'))&&(strcmp(mode.symmetry,'Cylindrical-X')))
yc = (y1+y2+y3)/3;
Se1=2*pi*yc.*Se1;
Se2=2*pi*yc.*Se2;
Se3=2*pi*yc.*Se3;
G1 = 2*pi*yc.*Gji(y2,y3,mode.tolGji); 
G2 = 2*pi*yc.*Gji(y3,y1,mode.tolGji); 
G3 = 2*pi*yc.*Gji(y1,y2,mode.tolGji); end
% =============================================================================================100
% here we define the logical masks for dirichlet boundary conditions;
% these masks are used when assembling the model equations:
%
%    1) maskr -> Poisson equation
%    2) masks -> elec/hole continuity equations, trap equations
%    3) maskt -> current equations
%
maskr=true(1,nn); masks=true(1,nn); maskt=false(1,nn); 
for ic=1:nc; ii=(contact==ic); % loop on contacts
switch geom.contact_type(ic)
case 1 % ohmic contact
maskr(ii)=0; % Poisson equation 
masks(ii)=0; % elec/hole equations
case 2 % Schottky contact
maskr(ii)=0; % Poisson equation 
case 3 % Schottky contact, vsurf=Inf
maskr(ii)=0; % Poisson equation 
masks(ii)=0; % elec/hole equations
otherwise, error('ls_solve -> contact type unknown!'), end, 
if(ic<=nm), maskt(ii)=1; end, end
%
% =============================================================================================100
% here we define the pointers to the matrix entries 
ijr=[in1 in2 in3];
iir=[in1 in1 in1 in2 in2 in2 in3 in3 in3];
jjr=[in1 in2 in3 in1 in2 in3 in1 in2 in3];
%
mask_iir=maskr(iir); 
mask_ijr=maskr(ijr); 
%
% the pointers IIT and JJT are necessary for the current equations;
in=zeros(1,nn);
for ic=1:nm; 
    ii=find(contact==1); % contacts
    in(ii)=rr+(ic-1)+1;
end
% 
iir=iir(mask_iir); jjr=jjr(mask_iir);
ijr=ijr(mask_ijr); 
%
% =============================================================================================100
% Loading material parameters and previous solution
% =============================================================================================100
% loading material parameters
dop_a=mesh.dop_a; dop_d=mesh.dop_d; epsxx=mesh.epsxx;
Nc=mesh.Nc; Nv=mesh.Nv; Eg=mesh.Eg;
%
% loading neutrality solution
phi_neutr=mode.phi_neutr; 
%
% loading previous solution
phi=uvet(1:nn); % normalized electric potential ( TO BE FIXED)
i_dd=0; % contact currents
v_dd=0; % contact voltages
uvet(qq+1)=v_dd;
% =============================================================================================100
% band diagram, eV
mode.ecb = NaN*zeros(1,nn); 
mode.evb = NaN*zeros(1,nn); 
mode.EFn = zeros(1,nn);
mode.EFp = zeros(1,nn);
mode.Efi = NaN*zeros(1,nn);
%
mode.ecb(iq) = - phi(iq) + Eg(iq)/2 - Vt(iq)./2.*log(Nv(iq)./Nc(iq)) + mesh.phi_r(iq);
mode.evb(iq) = mode.ecb(iq) - Eg(iq);
mode.Efi(iq) = (mode.ecb(iq) + mode.evb(iq))/2 + Vt(iq)./2.*log(Nv(iq)./Nc(iq));

%
if((isfield(mode,'stats'))&&(strcmp(mode.stats,'Fermi'))) % Fermi ---------
%
nF  = mesh.Nc.*ferdr((mode.EFn - mode.ecb)./Vt,1/2); % 1/cm^3;
dnF = mesh.Nc.*ferdr((mode.EFn - mode.ecb)./Vt,-1/2);
nB =  mesh.Nc.*exp  ((mode.EFn - mode.ecb)./Vt); % 1/cm^3; 
gamman = nF./nB; gamman(not(iq)) = 1;
dgamman = (1 - nF./dnF)./nB;
%
pF  = mesh.Nv.*ferdr((- mode.EFp + mode.evb)./Vt,1/2); % 1/cm^3;
dpF = mesh.Nv.*ferdr((- mode.EFp + mode.evb)./Vt,-1/2);
pB =  mesh.Nv.*exp  ((- mode.EFp + mode.evb)./Vt); % 1/cm^3; 
gammap = pF./pB; gammap(not(iq)) = 1;
dgammap = (1 - pF./dpF)./pB;
%
else % Boltzmann ----------------------------------------------------------
%    
nB =  mesh.Nc.*exp  ((mode.EFn - mode.ecb)./Vt); % 1/cm^3; 
nF  = nB; % 1/cm^3;
gamman = ones(size(nB));
dgamman = zeros(size(nB));
%
pB =  mesh.Nv.*exp  ((- mode.EFp + mode.evb)./Vt); % 1/cm^3; 
pF  = pB; % 1/cm^3;
gammap = ones(size(pB));
dgammap = zeros(size(pB));
end % ------------
%
% =============================================================================================100
phi=reshape(phi,1,nn);
% =============================================================================================100
% compute electric field (in time domain), V/cm
mode.efieldx=zeros(1,nt); mode.efieldy=zeros(1,nt);
mode.efieldx=-phi*mesh.gradx; 
mode.efieldy=-phi*mesh.grady; 
% =============================================================================================100
% assem Poisson equation 
% =============================================================================================100
% #############################################################################################100
if((isfield(mode,'ionization'))&&(strcmp(mode.ionization,'Incomplete')))
    gD = 2; n1 = mesh.Nc.*exp(-mesh.DeltaEd./Vt);
    dop_dp  =    dop_d./( 1 + gD.*nF./(n1.*gamman)); % ionized donor density, 1/cm^3
    ddop_dp =  - dop_d./((1 + gD.*nF./(n1.*gamman)).^2) ...
        .*gD.*((n1.*gamman - nF.*n1.*dgamman)./((n1.*gamman).^2));
    dop_dp(dop_d==0)=0;
    ddop_dp(dop_d==0)=0;
    %
    gA = 4; p1 = mesh.Nv.*exp(-mesh.DeltaEa./Vt);
    dop_am  =    dop_a./( 1 + gA.*pF./(p1.*gammap)); % ionized donor density, 1/cm^3
    ddop_am =  - dop_a./((1 + gA.*pF./(p1.*gammap)).^2) ...
        .*gA.*((p1.*gammap - pF.*p1.*dgammap)./((p1.*gammap).^2));
    dop_am(dop_a==0)=0;
    ddop_am(dop_a==0)=0;
    %
else dop_dp = dop_d;
     dop_am = dop_a; end
%
% save ionized dopant profiles 
mode.dop_dp = dop_dp;
mode.dop_am = dop_am;
% =============================================================================================100
it=ismember(triangle(4,:),find(geom.reg2contact)); % triangles inside conductors
M11 = ( s3./l3.*G3 + s2./l2.*G2).*epsxx; M11(it)=0;
M12 = (-s3./l3.*G3             ).*epsxx; M12(it)=0;
M13 = (-s2./l2.*G2             ).*epsxx; M13(it)=0;
M22 = ( s3./l3.*G3 + s1./l1.*G1).*epsxx; M22(it)=0;
M23 = (-s1./l1.*G1             ).*epsxx; M23(it)=0;
M33 = ( s1./l1.*G1 + s2./l2.*G2).*epsxx; M33(it)=0;
MM=[M11 M12 M13 M12 M22 M23 M13 M23 M33];
%
Kmat0=sparse(iir,jjr,MM(mask_iir),neq,neq);
% add current equations
Jmat0=sparse(neq,neq); 
% =============================================================================================100
qelNorm=qel*mode.CarrierNorm;
mode.qelNorm=qelNorm;
qelNorm2D=qel*mode.CarrierNorm2D;
mode.qelNorm2D=qelNorm2D;
%
% doping charge in Poisson rhs
MM=qelNorm.*[Se1.*dop_dp(in1) Se2.*dop_dp(in2) Se3.*dop_dp(in3)];
tvet=sparse(ijr,1,-MM(mask_ijr),neq,1); 
%
MM=qelNorm.*[Se1.*dop_am(in1) Se2.*dop_am(in2) Se3.*dop_am(in3)];
tvet=tvet+sparse(ijr,1,MM(mask_ijr),neq,1); 
%
% electron and hole charge in Poisson rhs (equilibrium)
if((isfield(mode,'stats'))&&(strcmp(mode.stats,'Fermi')))    
dnF =   mesh.Nc./Vt.*ferdr(( mode.EFn - mode.ecb)./Vt,-1/2);  % 1/cm^3;
dpF = - mesh.Nv./Vt.*ferdr((- mode.EFp + mode.evb)./Vt,-1/2); % 1/cm^3;
else dnF =   mesh.Nc./Vt.*exp((  mode.EFn - mode.ecb)./Vt);
     dpF = - mesh.Nv./Vt.*exp((- mode.EFp + mode.evb)./Vt); end 
%
nF(~iq)=0; dnF(~iq)=0; pF(~iq)=0; dpF(~iq)=0;
%
MM=qelNorm.*[Se1.*pF(in1) Se2.*pF(in2) Se3.*pF(in3)];
tvet=tvet+sparse(ijr,1,-MM(mask_ijr),neq,1); 
%
MM=qelNorm.*[Se1.*nF(in1) Se2.*nF(in2) Se3.*nF(in3)];
tvet=tvet+sparse(ijr,1, MM(mask_ijr),neq,1); 
%
MM=qelNorm.*[dpF(in1).*Se1 dpF(in2).*Se2 dpF(in3).*Se3];
Jmat0=Jmat0+sparse(ijr,ijr,-MM(mask_ijr),neq,neq);
%
MM=qelNorm.*[dnF(in1).*Se1 dnF(in2).*Se2 dnF(in3).*Se3];
Jmat0=Jmat0+sparse(ijr,ijr, MM(mask_ijr),neq,neq);
%
if((isfield(mode,'ionization'))&&(strcmp(mode.ionization,'Incomplete')))
MM=qelNorm.*[ddop_dp(in1).*dnF(in1).*Se1 ddop_dp(in2).*dnF(in2).*Se2 ddop_dp(in3).*dnF(in3).*Se3];
Jmat0=Jmat0+sparse(ijr,ijr,-MM(mask_ijr),neq,neq);
%
MM=qelNorm.*[ddop_am(in1).*dpF(in1).*Se1 ddop_am(in2).*dpF(in2).*Se2 ddop_am(in3).*dpF(in3).*Se3];
Jmat0=Jmat0+sparse(ijr,ijr,MM(mask_ijr),neq,neq); end
%
mode.elec(mesh.iq) = nF(mesh.iq);
mode.hole(mesh.iq) = pF(mesh.iq);
% =============================================================================================100
% Add circuit equations V-Z*I-V0=0
% =============================================================================================100
Zmat=sparse(mode.Zmat);
Jmat2=sparse(neq,neq);
Jmat2=Jmat2+sparse(qq+1,qq+1,1,neq,neq); % V
%
for ic=1:nm;
    ii=(ic-1)+1;
    for jc=1:nm;
        jj=(jc-1)+1;
        Jmat2(qq+ii,rr+jj)=Zmat(ii,jj);
    end
end % -Z*I
tvet(qq+(1:nm))=-v0_ls; % -V0
%
% =============================================================================================100
% Add current equations
% =============================================================================================100
% Using -qel instead of -1 because of equations normalization by qel
Kmat0=Kmat0+sparse(rr+(1:nm),rr+(1:nm),-1,neq,neq);
%
% =============================================================================================100
% Assembling quantum corrections and related jacobians
% =============================================================================================100
if(mode.qflg)
    N2Dtot=0; P2Dtot=0;
    for indQW=1:mesh.NMQW
        %
        % Basic parameters: geometry
        nnQW = mesh.nnxQW{indQW};
        inQW = mesh.inMQW{indQW}; % index to QW nodes
        %
        xQW = mesh.node(1,inQW);
        LeQW = xQW(2:end) - xQW(1:end-1); % edge length, cm
        Lp = zeros(1,nnQW); % box length, cm
        Lp(1:(nnQW-1)) = LeQW/2; Lp(2:nnQW) = Lp(2:nnQW) + LeQW/2;
        %
        % % debug: check Lp1,Lp2: integral of the constant function on the domain
        % %        this is expected to be equal to the length (z) of the QW
        % sum(diag(sparse([inQW1 inQW2],[inQW1 inQW2],[Lp1 Lp2],nn,nn)),1)
        %
        % Basic parameters: physics
        s_LoadConstants
        Vt2D = Vt(inQW);
        %
        DeltaEc = mesh.DeltaEcQW{indQW};
        DeltaEv = mesh.DeltaEvQW{indQW};
        lambda_C=mesh.lambda_C{indQW};
        lambda_V=mesh.lambda_V{indQW};
        Psi_C2=mesh.Psi_C2{indQW};
        Psi_V2=mesh.Psi_V2{indQW};
        %
        NBound_C=length(lambda_C); % number of conduction band bound states
        NBound_V=length(lambda_V); % number of valence band bound states
        Ec0 = ones(NBound_C,1)*mode.ecb(inQW);
        Ev0 = ones(NBound_V,1)*mode.evb(inQW);
        Em = Ec0 - DeltaEc + lambda_C*ones(1,nnQW); % eigenvalue referred to bottom of QW
        Ed = Ev0 + DeltaEv - lambda_V*ones(1,nnQW); % eigenvalue referred to top of QB
        Nc2D = mesh.Nc2D{indQW};%*ones(1,nnQW);
        Nv2D = mesh.Nv2D{indQW};%*ones(1,nnQW);
        MVt2D_C=ones(NBound_C,1)*Vt2D;
        MVt2D_V=ones(NBound_V,1)*Vt2D;
        %
        % @@@@@@@@@ Thermodynamic equilibrium solver @@@@@@@@@
        % add 2D charge to r.h.s. of Poisson equation
        EFn2D = zeros(NBound_C,1)*mode.EFn(inQW);
        EFp2D = zeros(NBound_V,1)*mode.EFp(inQW);
        xnum=(EFn2D-Em)./MVt2D_C; xden=(EFn2D-Ec0)./MVt2D_C;
        if(mode.iTappo==1 | mode.iTappo==3 | mode.iTappo==6)
            [n2Dm,dn2Dm] = ferdr2D(xnum,xden);
        else
            [n2Dm,dn2Dm] = ferdr2D_0(xnum,xden);
        end
        n2Dm=Nc2D.*n2Dm;
        dn2Dm=Nc2D.*dn2Dm./MVt2D_C;
        n2DmPois = zeros(size(Nc2D,1),length(mesh.xgrid)); 
        dn2DmPois = n2DmPois;
        n2DmPois(:,1:length(n2Dm))=n2Dm;
        dn2DmPois(:,1:length(dn2Dm))=dn2Dm;
        N2D = reshape(Psi_C2.'*n2DmPois,1,nn);
        N2Dtot = N2Dtot + N2D;
        dN2D = reshape(Psi_C2.'*dn2DmPois,1,nn);
        n2D = sum(n2Dm,1);
        mode.n2Di{indQW}=sum(n2Dm,1);
        %
        xnum=(Ed-EFp2D)./MVt2D_V; xden=(Ev0-EFp2D)./MVt2D_V;
        if(mode.iTappo==1 | mode.iTappo==3 | mode.iTappo==6)
            [p2Dm,dp2Dm] = ferdr2D(xnum,xden);
        else
            [p2Dm,dp2Dm] = ferdr2D_0(xnum,xden);
        end
        p2Dm=Nv2D.*p2Dm;
        dp2Dm=Nv2D.*dp2Dm./MVt2D_V;
        p2DmPois = zeros(size(Nv2D,1),length(mesh.xgrid)); 
        dp2DmPois = p2DmPois;
        p2DmPois(:,1:length(n2Dm))=p2Dm;
        dp2DmPois(:,1:length(dn2Dm))=dp2Dm;
        P2D = reshape(Psi_V2.'*p2DmPois,1,nn);
        P2Dtot = P2Dtot + P2D;
        dP2D = reshape(Psi_V2.'*(-dp2DmPois),1,nn);
        p2D = sum(p2Dm,1);
        mode.p2Di{indQW}=sum(p2Dm,1);
        %
        % assembling n2D charge in Poisson equation...
        MM = qelNorm2D.*[Se1.*N2D(in1) Se2.*N2D(in2) Se3.*N2D(in3)];
        tvet=tvet+sparse(ijr,1,MM(mask_ijr),neq,1);
        %
        % ...and jacobian.
        MM = qelNorm2D.*[Se1.*dN2D(in1) Se2.*dN2D(in2) Se3.*dN2D(in3)];
        Jmat0=Jmat0+sparse(ijr,ijr,MM(mask_ijr),neq,neq);
        %
        % assembling p2D charge in Poisson equation...
        MM = qelNorm2D.*[Se1.*P2D(in1) Se2.*P2D(in2) Se3.*P2D(in3)];
        tvet=tvet+sparse(ijr,1,-MM(mask_ijr),neq,1);
        %
        % ...and jacobian.
        MM = qelNorm2D.*[Se1.*dP2D(in1) Se2.*dP2D(in2) Se3.*dP2D(in3)];
        Jmat0=Jmat0+sparse(ijr,ijr,-MM(mask_ijr),neq,neq);
        %
        mode.Em{indQW}=Em; mode.n2D{indQW}=n2D; mode.EFn2D{indQW}=EFn2D;
        mode.Ed{indQW}=Ed; mode.p2D{indQW}=p2D; mode.EFp2D{indQW}=EFp2D;
        %
        if(mode.iTappo==6) % computing escape times as post-processing
            WQW = mesh.vWMQW{indQW}; %-- cm
            %
            meffn2D = mesh.meffnMQW{indQW};
            N2 = 4*pi/h^2*m0*qel*(sum((meffn2D*ones(1,mesh.nnxQW{indQW})).*(Ec0-Em),1)/10000);
            nBar = nF(inQW)*WQW;
            nWL = n2D;
            rn = nWL./(nBar.*(1-nWL/N2));
            mode.rn{indQW} = rn;
            %
            meffp2D = mesh.meffpMQW{indQW};
            P2 = 4*pi/h^2*m0*qel*(sum((meffp2D*ones(1,mesh.nnxQW{indQW})).*(Ed-Ev0),1)/10000);
            pBar = pF(inQW)*WQW;
            pWL = p2D;
            rp = pWL./(pBar.*(1-pWL/P2));
            mode.rp{indQW} = rp;
        end
    end
    mode.N2D=N2Dtot; mode.P2D=P2Dtot;
end
% =============================================================================================100
%
%
%
%
%
% =============================================================================================100
% Apply boundary conditions, Poisson equation: phi-V-phi_r=0
% =============================================================================================100
ii=find(not(maskr));
Kmat0=Kmat0+sparse(ii,ii,1,neq,neq); % phi
%
for ic=1:nm;
    ii=find(contact==ic);
    Kmat0=Kmat0+sparse(ii,qq+(ic-1)+1,-1,neq,neq);
end % - V
%
for ic=1:nc; ii=find(contact==ic); % loop on contacts
switch geom.contact_type(ic)
case 1 % ohmic contact
tvet(ii)=-phi_neutr(ii); % -phi_r
case 2 % Schottky contact
tvet(ii)=(geom.workfun(ic)+mesh.phi_rr); % -phi_r
case 3 % Schottky contact, vsurf=Inf
tvet(ii)=(geom.workfun(ic)+mesh.phi_rr); % -phi_r
otherwise, error('contact type unknown!'), end, end
% =============================================================================================100
% Compute residual
% =============================================================================================100
rvet = (Kmat0 + Jmat2) * uvet.' + tvet;
% =============================================================================================100
%
%
%
% *********************************************************************************************100
function [G] = Gji(rj,ri,threshold)
G=1; %-- Excluding calculation of this function
% G = (rj-ri)./log(rj./ri).*2./(rj+ri);
% indDiag = abs(rj-ri)<threshold; % diagonal (origin excluded)
% G(indDiag) = 1;
% indAxis_rj = abs(ri)<threshold; % rj axis (ri=0)
% G(indAxis_rj) = threshold;
% indAxis_ri = abs(rj)<threshold; % ri axis (rj=0)
% G(indAxis_ri) = threshold;
%-- Origin: no lim exists!
% indZer = abs(rj)<threshold & abs(ri)<threshold; 
% G(indZer) = threshold; % arbitrary choice!
%
% *********************************************************************************************100
function [x]=invferdr(y,tol)
%
% y=ferdr(x,1/2), so
% n = NC*ferdr(x,1/2) = NC*y
%
xini = log(y) + y.*(64+0.05524*y.*(64+sqrt(y))).^(-1/4);
x=xini;
res = Inf;
%
while res>tol
x0 = x; 
x = x0 - (ferdr(x0,1/2)-y)./ferdr(x0,-1/2);
res = max(abs(ferdr(x,1/2)-y)./y); end
%
% *********************************************************************************************100
function [h]=ferdr(x,k)
%
% compute Fermi integrals of order k=1/2,-1/2,3/2
%
p1=[-1.253314128820d+0, -1.723663557701d+0, -6.559045729258d-1, -6.342283197682d-2, -1.488383106116d-5];
q1=[+1.000000000000d+0, +2.191780925980d+0, +1.605812955406d+0, +4.443669527481d-1, +3.624232288112d-2];
p2=[-3.133285305570d-1, -4.161873852293d-1, -1.502208400588d-1, -1.339579375173d-2, -1.513350700138d-5];
q2=[+1.000000000000d+0, +1.872608675902d+0, +1.145204446578d+0, +2.570225587573d-1, +1.639902543568d-2];
p3=[-2.349963985406d-1, -2.927373637547d-1, -9.883097588738d-2, -8.251386379551d-3, -1.874384153223d-5];
q3=[+1.000000000000d+0, +1.608597109146d+0, +8.275289530880d-1, +1.522322382850d-1, +7.695120475064d-3];
p4=[+1.0738127694d+0, +5.6003303660d+0, +3.6882211270d+0, +1.1743392816d+0, +2.3641935527d-1];
q4=[+1.0000000000d+0, +4.6031840667d+0, +4.3075910674d-1, +4.2151132145d-1, +1.1832601601d-2];
p5=[+6.78176626660d-1, +6.33124017910d-1, +2.94479651772d-1, +8.01320711419d-2, +1.33918212940d-2];
q5=[+1.00000000000d+0, +1.43740400397d-1, +7.08662148450d-2, +2.34579494735d-3, -1.29449928835d-5];
p6=[+1.1530213402d+0, +1.0591558972d+0, +4.6898803095d-1, +1.1882908784d-1, +1.9438755787d-2];
q6=[+1.0000000000d+0, +3.7348953841d-2, +2.3248458137d-2, -1.3766770874d-3, +4.6466392781d-5];
p7=[-8.222559330d-1, -3.620369345d+1, -3.015385410d+3, -7.049871579d+4, -5.698145924d+4];
q7=[+1.000000000d+0, +3.935689841d+1, +3.568756266d+3, +4.181893625d+4, +3.385138907d+5];
p8=[+8.2244997626d-1, +2.0046303393d+1, +1.8268093446d+3, +1.2226530374d+4, +1.4040750092d+5];
q8=[+1.0000000000d+0, +2.3486207659d+1, +2.2013483743d+3, +1.1442673596d+4, +1.6584715900d+5];
p9=[+2.46740023684d+0, +2.19167582368d+2, +1.23829379075d+4, +2.20667724968d+5, +8.49442920034d+5];
q9=[+1.00000000000d+0, +8.91125140619d+1, +5.04575669667d+3, +9.09075946304d+4, +3.89960915641d+5];
c1=1.77245385090551603D0; c2=0.88622692545275801D0;
c3=1.32934038817913702D0; c4=0.66666666666666667D0;
c5=0.40000000000000000D0;
%
ii=(x<=1); jj=(x>1)&(x<=4); kk=not(ii|jj); h=[];
%
if(k==-1/2), % Fermi integral of order -1/2
y=exp(x(ii));
h(ii)=y.*(c1+y.*(p1(1)+y.*(p1(2)+y.*(p1(3)+y.*(p1(4)+y.*p1(5)))))./ ...
   (q1(1)+y.*(q1(2)+y.*(q1(3)+y.*(q1(4)+y.*q1(5))))));
h(jj)=(p4(1)+x(jj).*(p4(2)+x(jj).*(p4(3)+x(jj).*(p4(4)+x(jj).*p4(5)))))./ ...
   (q4(1)+x(jj).*(q4(2)+x(jj).*(q4(3)+x(jj).*(q4(4)+x(jj).*q4(5)))));
y=1./x(kk).^2;
h(kk)=sqrt(x(kk)).*(2+y.*(p7(1)+y.*(p7(2)+y.*(p7(3)+y.*(p7(4)+y.*p7(5)))))./ ...
   (q7(1)+y.*(q7(2)+y.*(q7(3)+y.*(q7(4)+y.*q7(5))))));
%
elseif(k==1/2), % Fermi integral of order 1/2
y=exp(x(ii));
h(ii)=y.*(c2+y.*(p2(1)+y.*(p2(2)+y.*(p2(3)+y.*(p2(4)+y.*p2(5)))))./ ...
   (q2(1)+y.*(q2(2)+y.*(q2(3)+y.*(q2(4)+y.*q2(5))))));
h(jj)=(p5(1)+x(jj).*(p5(2)+x(jj).*(p5(3)+x(jj).*(p5(4)+x(jj).*p5(5)))))./ ...
   (q5(1)+x(jj).*(q5(2)+x(jj).*(q5(3)+x(jj).*(q5(4)+x(jj).*q5(5)))));
y=1./x(kk).^2;
h(kk)=x(kk).*sqrt(x(kk)).*(c4+y.*(p8(1)+y.*(p8(2)+y.*(p8(3)+y.*(p8(4)+y.*p8(5)))))./ ...
   (q8(1)+y.*(q8(2)+y.*(q8(3)+y.*(q8(4)+y.*q8(5))))));
%
elseif(k==3/2), % Fermi integral of order 3/2
y=exp(x(ii));
h(ii)=y.*(c3+y.*(p3(1)+y.*(p3(2)+y.*(p3(3)+y.*(p3(4)+y.*p3(5)))))./ ...
   (q3(1)+y.*(q3(2)+y.*(q3(3)+y.*(q3(4)+y.*q3(5))))));
h(jj)=(p6(1)+x(jj).*(p6(2)+x(jj).*(p6(3)+x(jj).*(p6(4)+x(jj).*p6(5)))))./ ...
   (q6(1)+x(jj).*(q6(2)+x(jj).*(q6(3)+x(jj).*(q6(4)+x(jj).*q6(5)))));
y=1./x(kk).^2;
h(kk)=x(kk).^2.*sqrt(x(kk)).*(c5+y.*(p9(1)+y.*(p9(2)+y.*(p9(3)+y.*(p9(4)+y.*p9(5)))))./ ...
   (q9(1)+y.*(q9(2)+y.*(q9(3)+y.*(q9(4)+y.*q9(5))))));
end
%
h=1/gamma(k+1)*h;
% *********************************************************************************************100
