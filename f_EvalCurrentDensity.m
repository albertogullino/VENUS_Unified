
function [Jn_x,Jn_y,Jp_x,Jp_y] = f_EvalCurrentDensity(geom,mesh,mode)

if mode.Elementi==1
    [Jn_x,Jn_y,Jp_x,Jp_y] = f_EvalCurrentDensityElementi(geom,mesh,mode) ; 
    return
end
%
% function [Jn_x,Jn_y,Jp_x,Jp_y] = f_EvalCurrentDensity(geom,mesh,mode)
%
% The idea exploited to compute these currents is to treat each node as a
% contact. This can be obtained by assembling the contributions from a
% row/column (for y and x directed current components, respectively), but
% not the adjacent. In other words, only one row/column every two is
% assembled. This can be also seen as assembling half of the control box
% contribution, so that, instead of requiring the "null divergence" on a
% box, we require "half of it". And then normalize by the appropriate edge
% length.
%
% Any other approach as of this moment failed: all of them rely on
% quasi-Fermi levels (not supporting Fermi-Dirac distribution or
% interpolation of some quantity (carriers, mobilities, or whatever).
% Several discontinuities or spikes appear, from both these approaches.
%
% Alberto Tibaldi, 01/03/2017

% =============================================================================================100
% Loading and computing geometry parameters
% =============================================================================================100
%
s_LoadConstants
Vt=kB.*mesh.T/qel;
th_current=0;
%
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
%
xm1=1/2*(x2+x3); xm2=1/2*(x3+x1); xm3=1/2*(x1+x2);
ym1=1/2*(y2+y3); ym2=1/2*(y3+y1); ym3=1/2*(y1+y2);
% split points
xcc=1./2.*(x2.*f32xy+x1.*f32xy-y2.*f23yy-x3.*f23xy-x1.*f23xy+y3.*f23yy)./(f32xy-f23xy);
ycc=1./2.*(x2.*f23xx-y2.*f23xy-y1.*f23xy-x3.*f23xx+y3.*f32xy+y1.*f32xy)./(f32xy-f23xy);
% triangle area
Se=abs(r2x.*r1y-r2y.*r1x)./2;
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
% Total number of triangles: 2*(Nr_x*Nr_y)
Nr_x=mesh.nnx-1; % number of rectangles along x direction
Nr_y=mesh.nny-1; % number of rectangles along y direction
Nr=Nr_x*Nr_y;
%
% =============================================================================================100
% Loading solution database
% =============================================================================================100
nmodes=mode.nmodes; % number of optical modes (from VELM)
if(isfield(mesh,'nnxQW'))
nnQW=mesh.nnxQW;          % number of lateral points in the QW (x-directed)
else
nnQW=0;
end
nt=mesh.nt;               % number of triangles
nm=geom.nm;               % number of active contacts 
nc=geom.nc;               % number of contacts 
nl=mode.ntrap;            % number of trap levels 
%
nn=mesh.nn;               % -> elec eqs.
pp=nn+mode.nflg*nn;       % -> hole eqs.
tt=pp+mode.pflg*nn;       % -> trap eqs.
qq=tt+mode.tflg*nl*nn;    % -> circuit eqs.
rr=qq+nm;                 % -> current eqs.
%
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
masks(not(iq))=0; % continuity equations are eliminated outside semiconductors; 
maskt(not(iq))=0; % current equations are eliminated outside semiconductors; 
% =============================================================================================100
% here we define the pointers to the matrix entries 
ijr=[in1 in2 in3]; ijs=ijr;
iir=[in1 in1 in1 in2 in2 in2 in3 in3 in3]; iis=iir; iiu=iir;
jjr=[in1 in2 in3 in1 in2 in3 in1 in2 in3]; jjs=jjr; jjt=jjr;
%
mask_iir=maskr(iir); 
mask_iis=masks(iir); 
mask_iit=maskt(iir); 
%
mask_ijr=maskr(ijr); 
mask_ijs=masks(ijr); 
%
% the pointers IIT and JJT are necessary for the current equations;
in=zeros(1,nn);
for ic=1:nm; 
    ii=find(contact==1); % contacts
    in(ii)=rr+(ic-1)+1;
end
iit=in(iir);  
% 
iir=iir(mask_iir); jjr=jjr(mask_iir);
iis=iis(mask_iis); jjs=jjs(mask_iis);
iit=iit(mask_iit); jjt=jjt(mask_iit);
%
ijr=ijr(mask_ijr); 
ijs=ijs(mask_ijs);
%
dop_a=mesh.dop_a; dop_d=mesh.dop_d; epsxx=mesh.epsxx;
mobn0=mesh.mobn0_t; mobp0=mesh.mobp0_t; Nc=mesh.Nc; Nv=mesh.Nv; 
nint=mesh.nint; Eg=mesh.Eg; affinity=mesh.affinity;
%
elec=mode.elec; hole=mode.hole; phi=mode.phi;
phi1=phi(in1); phi2=phi(in2); phi3=phi(in3);
%
if((isfield(mode,'stats'))&&(strcmp(mode.stats,'Fermi'))) % Fermi ---------
%
if(mode.nflg), 
elec = abs(elec); % kludge: forcing electron charge to be positive    
tmp = elec./Nc; 
mode.EFn = mode.ecb + Vt.*f_invferdr(tmp,mode.tolinvferdr); end
nF  = mesh.Nc.*ferdr((mode.EFn - mode.ecb)./Vt,1/2); % 1/cm^3;
dnF = mesh.Nc.*ferdr((mode.EFn - mode.ecb)./Vt,-1/2);
nB =  mesh.Nc.*exp  ((mode.EFn - mode.ecb)./Vt); % 1/cm^3; 
gamman = nF./nB; gamman(not(iq)) = 1;
dgamman = (1 - nF./dnF)./nB;

if(not(mode.firstrun))
ii = tmp<1e-9; dgamman(ii) = - sqrt(2)/4./mesh.Nc(ii); % non-degenerate limit
end 
dlGn = dgamman./gamman; dlGn(not(iq)) = 0;
dlGn1 = dlGn(in1); dlGn2 = dlGn(in2); dlGn3 = dlGn(in3);
%
if(mode.pflg),    
hole = abs(hole); % kludge: forcing hole charge to be positive
tmp = hole./Nv; 
mode.EFp = mode.evb - Vt.*f_invferdr(tmp,mode.tolinvferdr); end
pF  = mesh.Nv.*ferdr((- mode.EFp + mode.evb)./Vt,1/2); % 1/cm^3;
dpF = mesh.Nv.*ferdr((- mode.EFp + mode.evb)./Vt,-1/2);
pB =  mesh.Nv.*exp  ((- mode.EFp + mode.evb)./Vt); % 1/cm^3; 
gammap = pF./pB; gammap(not(iq)) = 1;
dgammap = (1 - pF./dpF)./pB;
if(not(mode.firstrun))
ii = tmp<1e-9; dgammap(ii) = - sqrt(2)/4./mesh.Nv(ii); % non-degenerate limit
end  
dlGp = dgammap./gammap; dlGp(not(iq)) = 0;
dlGp1 = dlGp(in1); dlGp2 = dlGp(in2); dlGp3 = dlGp(in3); 
%
else % Boltzmann ----------------------------------------------------------
%    
if(mode.nflg), 
elec = abs(elec); % kludge: forcing electron charge to be positive    
tmp = elec./Nc; 
mode.EFn = mode.ecb + Vt.*log(tmp); end
nB =  mesh.Nc.*exp  ((mode.EFn - mode.ecb)./Vt); % 1/cm^3; 
nF  = nB; % 1/cm^3;
dnF = nF;
gamman = ones(size(nB));
dgamman = zeros(size(nB));
dlGn = dgamman;
dlGn1 = dlGn(in1); dlGn2 = dlGn(in2); dlGn3 = dlGn(in3);
%
if(mode.pflg),    
hole = abs(hole); % kludge: forcing hole charge to be positive
tmp = hole./Nv; 
mode.EFp = mode.evb - Vt.*log(tmp); end
pB =  mesh.Nv.*exp  ((- mode.EFp + mode.evb)./Vt); % 1/cm^3; 
pF  = pB; % 1/cm^3;
dpF = pF;
gammap = ones(size(pB));
dgammap = zeros(size(pB));
dlGp = 1./gammap.*dgammap; dlGp(not(iq)) = 0;
dlGp1 = dlGp(in1); dlGp2 = dlGp(in2); dlGp3 = dlGp(in3); end % ------------
%
offset_n=zeros(1,nn); offset_p=zeros(1,nn);
offset_n(iq)=affinity(iq)+Vt(iq).*log(Nc(iq));
offset_p(iq)=affinity(iq)-Vt(iq).*log(Nv(iq))+Eg(iq);
offset_n1=offset_n(in1); offset_n2=offset_n(in2); offset_n3=offset_n(in3);
offset_p1=offset_p(in1); offset_p2=offset_p(in2); offset_p3=offset_p(in3);
Vt1 = Vt(in1); Vt2 = Vt(in2); Vt3 = Vt(in3);
%
[Dn, Dp, dDndphi1, dDndphi2, dDndphi3, dDpdphi1, dDpdphi2, dDpdphi3] = ...
   mobility(mode,mobn0,mobp0,Se,r1x,r2x,r3x,r1y,r2y,r3y);
% the mobility function should be called for each material
% here we set to zero the mobility outside semiconductors
Dn(it)=0; Dp(it)=0;
% compute mobility from Einstein relation
mesh.mobn_t=Dn./mode.Vt_tr; mesh.mobp_t=Dp./mode.Vt_tr; % save field-dependent mobilities 
%
% =============================================================================================100
% Computing electron current densities
% =============================================================================================100
delta12 = (phi1 + offset_n1 + Vt1.*log(gamman(in1)))./Vt1 - (phi2 + offset_n2 + Vt2.*log(gamman(in2)))./Vt2;
delta23 = (phi2 + offset_n2 + Vt2.*log(gamman(in2)))./Vt2 - (phi3 + offset_n3 + Vt3.*log(gamman(in3)))./Vt3;
delta31 = (phi3 + offset_n3 + Vt3.*log(gamman(in3)))./Vt3 - (phi1 + offset_n1 + Vt1.*log(gamman(in1)))./Vt1;
%
B12=bern(delta12);
B21=bern(-delta12);
B23=bern(delta23); 
B32=bern(-delta23);
B31=bern(delta31);  
B13=bern(-delta31);
%
% =============================================================================================100
% Computation of Jn_y
% =============================================================================================100
% triangle side lengths
l1=sqrt(r1x.^2+r1y.^2); 
l2=sqrt(r2x.^2+r2y.^2); 
l3=sqrt(r3x.^2+r3y.^2);
%
% edges for flux computations
s1=sqrt((xm1-xcc).^2+(ym1-ycc).^2); 
s2=sqrt((xm2-xcc).^2+(ym2-ycc).^2); 
s3=sqrt((xm3-xcc).^2+(ym3-ycc).^2); 
%
% Finding the triangles to be assembled
it=zeros(1,mesh.nt);
ipoint=[1:Nr_x,Nr+1:Nr+Nr_x];
it(ipoint)=1;
for ind=2:2:Nr_y
    offset=(ind-1)*Nr_x;
    itemp=offset+ipoint;
    it(itemp)=1;
end
%
s1(not(it))=0; 
s2(not(it))=0; 
s3(not(it))=0;
%
l1=inf*ones(size(l1));
l2=inf*ones(size(l2));
%
Lp=zeros(1,mesh.nnx);
Lp(1:end-1)=diff(mesh.xgrid)/2;
Lp(2:end)=Lp(2:end)+diff(mesh.xgrid)/2;
LP=ones(mesh.nny,1)*Lp;
LP=reshape(LP,1,mesh.nn);
%
% Not assembling edge 1 (oblique): it does not contribute, since s1=0
M11=qel.*Dn.*(s3./l3.*B12+s2./l2.*B13);
M12=-qel.*Dn.*s3./l3.*B21;
M13=-qel.*Dn.*s2./l2.*B31;
M21=-qel.*Dn.*s3./l3.*B12;
M22=qel.*Dn.*(s3./l3.*B21);
M31=-qel.*Dn.*s2./l2.*B13;
M33=qel.*Dn.*(s2./l2.*B31);
M23=zeros(size(s1));
M32=zeros(size(s1));
MM=[M11 M12 M13 M21 M22 M23 M31 M32 M33];
Jnmat_y=sparse(iis,jjs,MM(mask_iis),nn,nn); % Assembly electron eqs.
Jn_y=(Jnmat_y*elec.').'./LP;
JN_Y=reshape(Jn_y,mesh.nny,mesh.nnx);
JN_Y(2:2:end,:)=-JN_Y(2:2:end,:); 
Jn_y=reshape(JN_Y,1,mesh.nn);
%
% =============================================================================================100
% Computation of Jn_x
% =============================================================================================100
% triangle side lengths
l1=sqrt(r1x.^2+r1y.^2); 
l2=sqrt(r2x.^2+r2y.^2); 
l3=sqrt(r3x.^2+r3y.^2);
%
% edges for flux computations
s1=sqrt((xm1-xcc).^2+(ym1-ycc).^2); 
s2=sqrt((xm2-xcc).^2+(ym2-ycc).^2); 
s3=sqrt((xm3-xcc).^2+(ym3-ycc).^2); 
%
it=zeros(1,mesh.nt);
ipoint=[1:2:Nr_x,Nr+[1:2:Nr_x]];
for ind=1:Nr_y
    offset=(ind-1)*Nr_x;
    itemp=offset+ipoint;
    it(itemp)=1;
end
%
s1(not(it))=0; 
s2(not(it))=0; 
s3(not(it))=0;
%
l1=inf*ones(size(l1));
l3=inf*ones(size(l3));
%
Lp=zeros(1,mesh.nny);
Lp(1:end-1)=diff(mesh.ygrid)/2;
Lp(2:end)=Lp(2:end)+diff(mesh.ygrid)/2;
%
% Further filtering a posteriori (if th_current~=0) the currents from too
% narrow edges
indSpikes=find(diff(mesh.ygrid)<th_current); 
indSpikes=indSpikes+1; % lavorando con diff, aggiungo 1
Lp(indSpikes)=inf;
LP=Lp.'*ones(1,mesh.nnx);
LP=reshape(LP,1,mesh.nn);
%
% Not assembling edge 1 (oblique): it does not contribute, since s1=0
M11=qel.*Dn.*(s3./l3.*B12+s2./l2.*B13);
M12=-qel.*Dn.*s3./l3.*B21;
M13=-qel.*Dn.*s2./l2.*B31;
M21=-qel.*Dn.*s3./l3.*B12;
M22=qel.*Dn.*(s3./l3.*B21);
M31=-qel.*Dn.*s2./l2.*B13;
M33=qel.*Dn.*(s2./l2.*B31);
M23=zeros(size(s1));
M32=zeros(size(s1));
MM=[M11 M12 M13 M21 M22 M23 M31 M32 M33];
Jnmat_x=sparse(iis,jjs,MM(mask_iis),nn,nn); % Assembly electron eqs.
Jn_x=(Jnmat_x*elec.').'./LP;
JN_X=reshape(Jn_x,mesh.nny,mesh.nnx);
JN_X(:,1:2:end)=-JN_X(:,1:2:end); 
Jn_x=reshape(JN_X,1,mesh.nn);
%
% =============================================================================================100
% Computing hole current densities
% =============================================================================================100
delta12 = (phi1 + offset_p1 - Vt1.*log(gammap(in1)))./Vt1 - (phi2 + offset_p2 - Vt2.*log(gammap(in2)))./Vt2;
delta23 = (phi2 + offset_p2 - Vt2.*log(gammap(in2)))./Vt2 - (phi3 + offset_p3 - Vt3.*log(gammap(in3)))./Vt3;
delta31 = (phi3 + offset_p3 - Vt3.*log(gammap(in3)))./Vt3 - (phi1 + offset_p1 - Vt1.*log(gammap(in1)))./Vt1;
%
B12=bern(delta12);
B21=bern(-delta12);
B23=bern(delta23);
B32=bern(-delta23);
B31=bern(delta31);
B13=bern(-delta31);
%
% =============================================================================================100
% Computation of Jp_y
% =============================================================================================100
% triangle side lengths
l1=sqrt(r1x.^2+r1y.^2); 
l2=sqrt(r2x.^2+r2y.^2); 
l3=sqrt(r3x.^2+r3y.^2);
%
% edges for flux computations
s1=sqrt((xm1-xcc).^2+(ym1-ycc).^2); 
s2=sqrt((xm2-xcc).^2+(ym2-ycc).^2); 
s3=sqrt((xm3-xcc).^2+(ym3-ycc).^2); 
%
% Finding the triangles to be assembled
it=zeros(1,mesh.nt);
ipoint=[1:Nr_x,Nr+1:Nr+Nr_x];
it(ipoint)=1;
for ind=2:2:Nr_y
    offset=(ind-1)*Nr_x;
    itemp=offset+ipoint;
    it(itemp)=1;
end
%
s1(not(it))=0; 
s2(not(it))=0; 
s3(not(it))=0;
%
l1=inf*ones(size(l1));
l2=inf*ones(size(l2));
%
Lp=zeros(1,mesh.nnx);
Lp(1:end-1)=diff(mesh.xgrid)/2;
Lp(2:end)=Lp(2:end)+diff(mesh.xgrid)/2;
LP=ones(mesh.nny,1)*Lp;
LP=reshape(LP,1,mesh.nn);
%
s1(not(it))=0; 
s2(not(it))=0; 
s3(not(it))=0;
%
l1=inf*ones(size(l1));
l2=inf*ones(size(l2));
%
M11=qel.*Dp.*(s3./l3.*B21+s2./l2.*B31);
M12=-qel.*Dp.*s3./l3.*B12;
M13=-qel.*Dp.*s2./l2.*B13;
M21=-qel.*Dp.*s3./l3.*B21;
M22=qel.*Dp.*(s3./l3.*B12);
M31=-qel.*Dp.*s2./l2.*B31;
M33=qel.*Dp.*(s2./l2.*B13);
M23=zeros(size(s1));
M32=zeros(size(s1));
MM=[M11 M12 M13 M21 M22 M23 M31 M32 M33];
Jpmat_y=sparse(iis,jjs,MM(mask_iis),nn,nn);
Jp_y=(Jpmat_y*hole.').'./LP;
JP_Y=reshape(Jp_y,mesh.nny,mesh.nnx);
JP_Y(1:2:end,:)=-JP_Y(1:2:end,:); 
Jp_y=reshape(JP_Y,1,mesh.nn);
%
% =============================================================================================100
% Computation of Jp_x
% =============================================================================================100
% triangle side lengths
l1=sqrt(r1x.^2+r1y.^2); 
l2=sqrt(r2x.^2+r2y.^2); 
l3=sqrt(r3x.^2+r3y.^2);
%
% edges for flux computations
s1=sqrt((xm1-xcc).^2+(ym1-ycc).^2); 
s2=sqrt((xm2-xcc).^2+(ym2-ycc).^2); 
s3=sqrt((xm3-xcc).^2+(ym3-ycc).^2); 
%
it=zeros(1,mesh.nt);
ipoint=[1:2:Nr_x,Nr+[1:2:Nr_x]];
for ind=1:Nr_y
    offset=(ind-1)*Nr_x;
    itemp=offset+ipoint;
    it(itemp)=1;
end
%
s1(not(it))=0; 
s2(not(it))=0; 
s3(not(it))=0;
%
l1=inf*ones(size(l1));
l3=inf*ones(size(l3));
%
Lp=zeros(1,mesh.nny);
Lp(1:end-1)=diff(mesh.ygrid)/2;
Lp(2:end)=Lp(2:end)+diff(mesh.ygrid)/2;
%
% Further filtering a posteriori (if th_current~=0) the currents from too
% narrow edges
indSpikes=find(diff(mesh.ygrid)<th_current); 
indSpikes=indSpikes+1; % lavorando con diff, aggiungo 1
Lp(indSpikes)=inf;
LP=Lp.'*ones(1,mesh.nnx);
LP=reshape(LP,1,mesh.nn);
%
M11=qel.*Dp.*(s3./l3.*B21+s2./l2.*B31);
M12=-qel.*Dp.*s3./l3.*B12;
M13=-qel.*Dp.*s2./l2.*B13;
M21=-qel.*Dp.*s3./l3.*B21;
M22=qel.*Dp.*(s3./l3.*B12);
M31=-qel.*Dp.*s2./l2.*B31;
M33=qel.*Dp.*(s2./l2.*B13);
M23=zeros(size(s1));
M32=zeros(size(s1));
MM=[M11 M12 M13 M21 M22 M23 M31 M32 M33];
Jpmat_x=sparse(iis,jjs,MM(mask_iis),nn,nn);
Jp_x=(Jpmat_x*hole.').'./LP;
JP_X=reshape(Jp_x,mesh.nny,mesh.nnx);
JP_X(:,2:2:end)=-JP_X(:,2:2:end); 
JP_x=reshape(JP_X,1,mesh.nn);
Jp_x=reshape(JP_X,1,mesh.nn);

return




JY=JP_Y+JN_Y;
for indz=1:mesh.nny
    I(indz)=trapz(mesh.xgrid,2*pi*mesh.xgrid.*JY(indz,:))
end

figure,plot(mesh.ygrid,I)






% figure,plot3(mesh.node(1,:),mesh.node(2,:),(Jn_x),'*')%,set(gca,'zscale','log'),zlim([1e-15,1])

figure,surf(mesh.X,mesh.Y,(JP_Y+JN_Y)),xlabel('\rho'),ylabel('z'),title('J_z'),shading interp
figure,surf(mesh.X,mesh.Y,(JP_X+JN_X)),xlabel('\rho'),ylabel('z'),title('J_{\rho}'),shading interp


% return

XX=mesh.X(1:passo_y:end,1:passo_x:end);
YY=mesh.Y(1:passo_y:end,1:passo_x:end);
JX=JP_X(1:passo_y:end,1:passo_x:end)+JN_X(1:passo_y:end,1:passo_x:end);
JY=JP_Y(1:passo_y:end,1:passo_x:end)+JN_Y(1:passo_y:end,1:passo_x:end);

% figure,quiver(mesh.X,mesh.Y,JN_X,JN_Y,3)
figure
hold on
grid on
box on
quiver(XX,YY,JX,JY,arrow_normalization)
axis equal
