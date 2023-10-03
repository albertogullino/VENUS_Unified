function [Kmat0,Jmat0,Jmat2,uvet,rvet,mode,tvet,tvet3]=assem_GBT_Elementi(geom,mesh,mode,uvet,v0_dd)
%
% define constants
%
s_LoadConstants
Vt=mode.Vt; % electrical equivalent in temperature, V
% =============================================================================================100
% Loading solution database
% =============================================================================================100
nmodes=mode.nmodes; % number of optical modes (from VELM)
if(isfield(mesh,'nnxQW'))
    nnQW=mesh.nnxQW{1};          % number of lateral points in the QW (x-directed)
else
    nnQW=0;
end
NQW=mesh.NMQW;                          % number of quantum wells
nt=mesh.nt;                             % number of triangles
nm=geom.nm;                             % number of active contacts
nc=geom.nc;                             % number of contacts
nl=mode.ntrap;                          % number of trap levels
%
if mode.EqPermutazione==0
    %
    nn=mesh.nn;                   % -> elec eqs.
    pp=nn+mode.nflg*nn;           % -> hole eqs.
    tt=pp+mode.pflg*nn;           % -> trap eqs.
    qq=tt+mode.tflg*nl*nn;        % -> circuit eqs.
    rr=qq+nm;                     % -> current eqs.
    vv=rr+nm;                     % -> n2D eqs.
    ww=vv+length(mesh.xgrid)*NQW;               % -> p2D eqs.
    ss=ww+length(mesh.xgrid)*NQW;               % -> optical rate eqs.
else
    nn=mesh.nn;                             % -> elec eqs.
    pp=nn+mode.nflg*nn;                     % -> hole eqs.
    tt=pp+mode.pflg*nn;                     % -> trap eqs.
    vv=tt+mode.tflg*nl*nn;                  % -> n2D eqs
    ww=vv+length(mesh.xgrid)*NQW;           % -> p2D eqs
    ss=ww+length(mesh.xgrid)*NQW;           % -> optical rate eqs.
    qq=ss+nmodes;
    rr=qq+nm;
end
%
% total number of equations/unknowns
neq=nn+mode.nflg*nn+mode.pflg*nn+mode.tflg*nl*nn+2*nm+mode.oflg*not(mode.firstrun)*2*NQW*length(mesh.xgrid)+mode.oflg*nmodes;
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
if(mode.AugerExpulsion3D==1)
    SeCAV1=zeros(1,nt);
    SeCAV2=SeCAV1;
    SeCAV3=SeCAV1;
    SeCAV1(mesh.ITrCAV)=Se1(mesh.ITrCAV);
    SeCAV2(mesh.ITrCAV)=Se2(mesh.ITrCAV);
    SeCAV3(mesh.ITrCAV)=Se3(mesh.ITrCAV);
end
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
    if(mode.AugerExpulsion3D==1)
        SeCAV1=2*pi*xc.*SeCAV1;
        SeCAV2=2*pi*xc.*SeCAV2;
        SeCAV3=2*pi*xc.*SeCAV3;
    end
    G1 = 2*pi*xc.*Gji(x2,x3,mode.tolGji);
    G2 = 2*pi*xc.*Gji(x3,x1,mode.tolGji);
    G3 = 2*pi*xc.*Gji(x1,x2,mode.tolGji);
end
%
if((isfield(mode,'symmetry'))&&(strcmp(mode.symmetry,'Cylindrical-X')))
    yc = (y1+y2+y3)/3;
    Se1=2*pi*yc.*Se1;
    Se2=2*pi*yc.*Se2;
    Se3=2*pi*yc.*Se3;
    if(mode.AugerExpulsion3D==1)
        SeCAV1=2*pi*yc.*SeCAV1;
        SeCAV2=2*pi*yc.*SeCAV2;
        SeCAV3=2*pi*yc.*SeCAV3;
    end
    G1 = 2*pi*yc.*Gji(y2,y3,mode.tolGji);
    G2 = 2*pi*yc.*Gji(y3,y1,mode.tolGji);
    G3 = 2*pi*yc.*Gji(y1,y2,mode.tolGji);
end
% =============================================================================================100
% here we define the logical masks for dirichlet boundary conditions;
% these masks are used when assembling the model equations:
%
%    1) maskr -> Poisson equation
%    2) masks -> elec/hole continuity equations, trap equations
%    3) maskt -> current equations
%
maskr=true(1,nn); masks=true(1,nn); maskt=false(1,nn);
for ic=1:nc
    ii=(contact==ic); % loop on contacts
    switch geom.contact_type(ic)
        case 1 % ohmic contact
            maskr(ii)=0; % Poisson equation
            masks(ii)=0; % elec/hole equations
        case 2 % Schottky contact
            maskr(ii)=0; % Poisson equation
        case 3 % Schottky contact, vsurf=Inf
            maskr(ii)=0; % Poisson equation
            masks(ii)=0; % elec/hole equations
        otherwise, error('ls_solve -> contact type unknown!'),
    end
    if(ic<=nm)
        maskt(ii)=1;
    end
end
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
for ic=1:nm
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
% =============================================================================================100
% Loading material parameters and previous solution
% =============================================================================================100
% loading material parameters
if isfield(mode,'Fat_Dop')
    FatDop=mode.Fat_Dop;
else
    FatDop=1;
end
dop_a=FatDop*mesh.dop_a; dop_d=FatDop*mesh.dop_d; epsxx=mesh.epsxx;
mobn0=mesh.mobn0_t; mobp0=mesh.mobp0_t; Nc=mesh.Nc; Nv=mesh.Nv;
nint=mesh.nint;
nint_equil=mesh.nint_equil; %

Eg=mesh.Eg;
affinity=mesh.affinity;

if isfield(mode,'ireno')==1
    ireno=mode.ireno;
    if ireno==1
        Reno=2.7*10^(-8)*mode.dop_am.^(1/3);
        Eg=Eg-Reno;
        affinity=affinity-Reno*.4;
        %
        'EG',keyboard
    end
end
% loading neutrality solution
phi_neutr=mode.phi_neutr;
%
% loading previous solution
phi=uvet(1:nn); % electric potential
if mode.IdriveON==1
    i0_dd=v0_dd;    % current driving
    Ymat=mode.Zmat;
    
    v_ls=uvet(qq+1);  % contact voltages
    %   i_ls=-v_ls*(Ymat.')+i0_dd;  % contact currents
    %   uvet(rr+1)=i_ls;
else
    i_ls=uvet(rr+1); % contact currents
    v_ls=-i_ls*(mode.Zmat.')+v0_dd; % contact voltages
    %     uvet(qq+1)=v_ls;
end
if(mode.nflg)
    elec=abs(uvet(nn+(1:nn)));
    
end % electron concentration
if(mode.pflg)
    hole=abs(uvet(pp+(1:nn)));
    
end % hole concentration

if(mode.oflg)
    Pst=abs(uvet(ss+(1:nmodes)));
    
    for indQW=1:NQW
        n2Dc{indQW} = abs(uvet(vv+(indQW-1)*length(mesh.xgrid)+(1:nnQW)));
        p2Dc{indQW} = abs(uvet(ww+(indQW-1)*length(mesh.xgrid)+(1:nnQW)));
    end
end
% =============================================================================================100
% these vectors are necessary for the treatment of heterointerfaces

Vt1 = Vt(in1); Vt2 = Vt(in2); Vt3 = Vt(in3);
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
if(mode.nflg), mode.EFn(iq) = mode.ecb(iq) + Vt(iq).*log(elec(iq)./Nc(iq)); end
if(mode.pflg), mode.EFp(iq) = mode.evb(iq) - Vt(iq).*log(hole(iq)./Nv(iq)); end
%
if((isfield(mode,'stats'))&&(strcmp(mode.stats,'Fermi'))) % Fermi ---------
    %
    if(mode.nflg)
        elec = abs(elec); % kludge: forcing electron charge to be positive
        tmp = elec./Nc;
        mode.EFn = mode.ecb + Vt.*invferdr(tmp,mode.tolinvferdr);
    end
    nF  = Nc.*ferdr((mode.EFn - mode.ecb)./Vt,1/2); % 1/cm^3;
    dnF = Nc.*ferdr((mode.EFn - mode.ecb)./Vt,-1/2);
    nB =  Nc.*exp  ((mode.EFn - mode.ecb)./Vt); % 1/cm^3;
    gamman = nF./nB; gamman(not(iq)) = 1;
    dgamman = (1 - nF./dnF)./nB;
    if(not(mode.firstrun))
        ii = tmp<1e-9; dgamman(ii) = - sqrt(2)/4./Nc(ii); % non-degenerate limit
    end
    dlGn = dgamman./gamman; dlGn(not(iq)) = 0;
    dlGn1 = dlGn(in1); dlGn2 = dlGn(in2); dlGn3 = dlGn(in3);
    %
    if(mode.pflg)
        hole = abs(hole); % kludge: forcing hole charge to be positive
        tmp = hole./Nv;
        mode.EFp = mode.evb - Vt.*invferdr(tmp,mode.tolinvferdr);
    end
    pF  = Nv.*ferdr((- mode.EFp + mode.evb)./Vt,1/2); % 1/cm^3;
    dpF = Nv.*ferdr((- mode.EFp + mode.evb)./Vt,-1/2);
    pB =  Nv.*exp  ((- mode.EFp + mode.evb)./Vt); % 1/cm^3;
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
    if(mode.nflg)
        elec = abs(elec); % kludge: forcing electron charge to be positive
        tmp = elec./Nc;
        mode.EFn = mode.ecb + Vt.*log(tmp);
    end
    nB =  Nc.*exp((mode.EFn - mode.ecb)./Vt); % 1/cm^3;
    nF  = nB; % 1/cm^3;
    dnF = nF;
    gamman = ones(size(nB));
    dgamman = zeros(size(nB));
    dlGn = dgamman;
    dlGn1 = dlGn(in1); dlGn2 = dlGn(in2); dlGn3 = dlGn(in3);
    %
    if(mode.pflg)
        hole = abs(hole); % kludge: forcing hole charge to be positive
        tmp = hole./Nv;
        mode.EFp = mode.evb - Vt.*log(tmp);
    end
    pB =  Nv.*exp((- mode.EFp + mode.evb)./Vt); % 1/cm^3;
    pF  = pB; % 1/cm^3;
    dpF = pF;
    gammap = ones(size(pB));
    dgammap = zeros(size(pB));
    dlGp = 1./gammap.*dgammap; dlGp(not(iq)) = 0;
    dlGp1 = dlGp(in1); dlGp2 = dlGp(in2); dlGp3 = dlGp(in3);
end % ------------
%
% =============================================================================================100
phi=reshape(phi,1,nn);
phi1=phi(in1); phi2=phi(in2); phi3=phi(in3);
if(mode.nflg)
    elec=reshape(elec,1,nn);
    elec1=elec(in1); elec2=elec(in2); elec3=elec(in3);
end
if(mode.pflg)
    hole=reshape(hole,1,nn);
    hole1=hole(in1); hole2=hole(in2); hole3=hole(in3);
end
% =============================================================================================100
% Mobility models
[Dn, Dp, dDndphi1, dDndphi2, dDndphi3, dDpdphi1, dDpdphi2, dDpdphi3] = ...
    mobility(mode,mobn0,mobp0,Se,r1x,r2x,r3x,r1y,r2y,r3y);
% the mobility function should be called for each material
% here we set to zero the mobility outside semiconductors
Dn(it)=0; Dp(it)=0;
% compute mobility from Einstein relation
mode.mobn=Dn./mode.Vt_tr; mode.mobp=Dp./mode.Vt_tr; % save field-dependent mobilities
% =============================================================================================100
% GR models
% =============================================================================================100
R=zeros(1,nn); dRn=zeros(1,nn); dRp=zeros(1,nn);
RSRH=zeros(1,nn);
Rrad=zeros(1,nn);   PspBulk=0;
RAuger=zeros(1,nn);
GBTBT=zeros(1,nn);

RAuger_n3D = zeros(1,nn); RAuger_p3D = zeros(1,nn);
dRAuger_n3D_n3D = zeros(1,nn); dRAuger_n3D_p3D = zeros(1,nn);
dRAuger_p3D_n3D = zeros(1,nn); dRAuger_p3D_p3D = zeros(1,nn);

%
% add current equations
Jmat0=sparse(neq,neq);
tvet=sparse(neq,1);

if(not(mode.firstrun)) % At equilibrium, R=0
    for indGR=1:length(mode.GR)
        if(not(mode.nflg && mode.pflg))
            error('ls_solve -> select both elec and hole eqs. with GR models')
        end
        n  =  elec(iq); p  =  hole(iq);
        np = (n.*p-gamman(iq).*gammap(iq).*nint_equil(iq).^2);
        dnp_n = (p-dgamman(iq).*gammap(iq).*nint_equil(iq).^2);
        dnp_p = (n-dgammap(iq).*gamman(iq).*nint_equil(iq).^2);
        
        if strcmp(mode.GR{indGR},'SRH')==1
            % SRH recombination model
            n1 = gamman(iq).*nint_equil(iq).*exp(mesh.Etrap(iq)./Vt(iq));
            p1 = gammap(iq).*nint_equil(iq).*exp(-mesh.Etrap(iq)./Vt(iq));
            denSRH = mesh.taup(iq).*(n+n1)+mesh.taun(iq).*(p+p1);
            RSRH(iq) = np./denSRH;
            R(iq) = R(iq) + RSRH(iq);
            dRn(iq) = dRn(iq) + ((p-dgamman(iq).*gammap(iq).*nint_equil(iq).^2).*denSRH - ...
                mesh.taup(iq).*(1+dgamman(iq).*n1).*np)./(denSRH.^2);
            dRp(iq) = dRp(iq) + ((n-dgammap(iq).*gamman(iq).*nint_equil(iq).^2).*denSRH - ...
                mesh.taun(iq).*(1+dgammap(iq).*p1).*np)./(denSRH.^2);
        elseif strcmp(mode.GR{indGR},'rad')==1
            % radiative recombination model
            Rrad(iq) = mesh.brad(iq).*np;
            R(iq) = R(iq) + Rrad(iq);
            dRn(iq)=dRn(iq)+mesh.brad(iq).*dnp_n;
            dRp(iq)=dRp(iq)+mesh.brad(iq).*dnp_p;
        elseif strcmp(mode.GR{indGR},'Auger')==1
            % Auger recombination model
            %
            RAuger_n3D(iq) = mesh.Cnnp(iq).*n.*np;
            RAuger_p3D(iq) = mesh.Cppn(iq).*p.*np;
            
            dRAuger_n3D_n3D(iq) = 2*mesh.Cnnp(iq).*np;
            dRAuger_n3D_p3D(iq) = mesh.Cnnp(iq).*dnp_p;
            
            dRAuger_p3D_n3D(iq) = mesh.Cppn(iq).*dnp_n;
            dRAuger_p3D_p3D(iq) = 2*mesh.Cppn(iq).*np;
            %
            RAuger(iq) = RAuger_n3D(iq) + RAuger_p3D(iq);
            %
            R(iq) = R(iq) + RAuger(iq);
            dRn(iq) = dRn(iq) + dRAuger_n3D_n3D(iq) + dRAuger_p3D_n3D(iq);
            dRp(iq) = dRp(iq) + dRAuger_n3D_p3D(iq) + dRAuger_p3D_p3D(iq);
        elseif (strcmp(mode.GR{indGR},'BTBT')==1)
            %
            % Load the NEGF LUT
            if (abs(v0_dd)==0)
                LUT_GBTBT
            end
            %
            % BTBT generation rate, from NEGF
            %
            for indexBTJ = 1:mode.nBTJ
                %
                GBTBT=zeros(1,nn);dRnTJ=zeros(1,mesh.nn); dRpTJ=dRnTJ; dRphin=dRnTJ; dRphip=dRpTJ;
                %
                iTJfake=indexBTJ;
                %
                IBTJ=sort(unique([mesh.IBTJ{indexBTJ,:}]));
                ITJp_vet=mesh.iLBTJ{iTJfake,end};                     % right BTJ node, (no passivation)
                ITJn_vet=mesh.iRBTJ{iTJfake,1};                     % left BTJ node
                iTJp=ITJp_vet(1);
                iTJn=ITJn_vet(1);
                %
                iRBTJP=ITJn_vet;                                % RIGHT BTJ nodes, all the columns
                iLBTJP=ITJp_vet;                                % LEFT BTJ nodes, all the columns
                %
                IIBTJP=repmat(1:mesh.nnx,mesh.nny,1);   IIBTJP=IIBTJP(:).';
                %
                IRBTJP=iRBTJP(IIBTJP);  % needed to assemble the GBTBT derivatives only in the right and left nodes of the BTJ (as the QW: INMQWP)
                ILBTJP=iLBTJP(IIBTJP);
                %
                SeTJ1=zeros(1,nt); SeTJ2=SeTJ1; SeTJ3=SeTJ1;
                SeTJ1(unique([mesh.ITrBTJ{indexBTJ,:}]))=Se1(unique([mesh.ITrBTJ{indexBTJ,:}]));
                SeTJ2(unique([mesh.ITrBTJ{indexBTJ,:}]))=Se2(unique([mesh.ITrBTJ{indexBTJ,:}]));
                SeTJ3(unique([mesh.ITrBTJ{indexBTJ,:}]))=Se3(unique([mesh.ITrBTJ{indexBTJ,:}]));
                %
                % Rule of thumb:
                % - BTJ is pn: n --> mesh.inTJ(end); p --> mesh.inTJ(1)
                % - BTJ is np (2020Tibaldi_PRAPP): n --> mesh.inTJ(1); p --> mesh.inTJ(end)
                %
                Vint=mode.EFp(ITJp_vet)-mode.EFn(ITJn_vet);
                Vint(isnan(Vint))=0;
                mode.VTJ(indexBTJ,mode.ind_v0)=Vint(1);
                % Initialization of Vinterp derivatives for Jacobian matrix
                % dn -> right (end); dp -> left (1); 1 -> right; 2 -> left
                % Vt,ecb, evb, EFn, EFp should ALL be taken from "real"
                % nodes, and not only from the central section!
                % iTJend -> IBTJright; iTJ1 -> IBTJleft
                dVint_dn = - Vt(ITJn_vet)./Nc(ITJn_vet)./ferdr((mode.EFn(ITJn_vet)-mode.ecb(ITJn_vet))./Vt(ITJn_vet),-1/2);
                dVint_dp = - Vt(ITJp_vet)./Nv(ITJp_vet)./ferdr(-(mode.EFp(ITJp_vet)-mode.evb(ITJp_vet))./Vt(ITJp_vet),-1/2);
                dVint_dn(isnan(dVint_dn))=0; dVint_dp(isnan(dVint_dp))=0;
                % La derivata rispetto al potenziale si calcola usando EFn = Ec
                % + altro,
                % dphi1 --> dRphi1, TJ(end); dphin --> dRphi2, TJ(1)
                dVint_dphip = -1; % derivata rispetto a EFp (-1 per colpa di Ev = -q*phi + puffi)
                dVint_dphin = +1; % derivata rispetto a EFn (-1 per colpa di Ev = -q*phi + puffi, altro -1 perché ha segno -)
                
                % TIBALDI JACOBIAN GBTBT
                % lunghezza della regione dove spalmo il tasso
                LTJ = mesh.ygrid(iTJn)-mesh.ygrid(iTJp);
                
                pDeg=size(mode.cBTBT,2)-1;
                
                VV=zeros(length(ITJp_vet),pDeg+1);
                dVV=zeros(length(ITJp_vet),pDeg);
                %
                [~,iRagTJ]=min(abs(mesh.xgrid*1e4-mode.rAperture));
                %
                %                     for iRad=1:length(ITJp_vet)
                for iRad=1:iRagTJ
                    
                    [~,iT]=min(abs(mode.T_NEGF-mean(mesh.T(ITJp_vet(iRad):ITJn_vet(iRad)))));
                    cT=mode.cBTBT(iT,:); % fitting coefficient at a given temperature
                    %             cT=c;       % fitting of NEGF w/out temperature
                    %
                    % coefficienti della derivata del fit della corrente (3 values, polyfit 3)
                    dc =cT(1:pDeg).*[pDeg:-1:1];
                    
                    VV(iRad,:) = Vint(iRad).^[pDeg:-1:0];
                    dVV(iRad,:) = Vint(iRad).^[pDeg-1:-1:0];
                    
                end
                
                I_interp=10.^(VV*cT')-mode.I0_NEGF(iT); % A/cm^2 !!!!
                II=1/LTJ*(-I_interp/qel);
                II(iRagTJ+1:end)=0;
                GBTBT(IBTJ)=II(IIBTJP(IBTJ))/mode.CarrierNorm;
                %   'BTBT',keyboard
                ITJ=I_interp*mode.AreaOx;   % A
                RTJ=abs(Vint)'./abs(ITJ);   % Ohm
                mode.sigmaTJ(indexBTJ)=1./RTJ(1)./LTJ;   % S/cm (compatible with f_EvalHeatingsTerms.m)
                
                JTJ=I_interp;
                mode.HeatTJ=abs(JTJ.*Vint'./repmat(LTJ,mesh.nnx,1));
                
                dI_interp=log(10)*10.^(VV*cT').*(dVV*dc');
                dII=1/LTJ*(-dI_interp/qel)'/mode.CarrierNorm;
                
                dRnTJ(IBTJ) = dII(IIBTJP(IBTJ)).*dVint_dn(IIBTJP(IBTJ));
                dRpTJ(IBTJ) = dII(IIBTJP(IBTJ)).*dVint_dp(IIBTJP(IBTJ));
                dRphin(IBTJ) = dII(IIBTJP(IBTJ)).*dVint_dphin;
                dRphip(IBTJ) = dII(IIBTJP(IBTJ)).*dVint_dphip;
                %
                % GBTBT derivatives assembled
                %
                % Assembling in electron continuity equation
                % Recombination term
                MM = qel.*[SeTJ1.*GBTBT(in1) SeTJ2.*GBTBT(in2) SeTJ3.*GBTBT(in3)];
                tvet = tvet + sparse(ijs+nn,1,MM(mask_ijs),neq,1);
                %
                % Derivative of recombination term w.r.t. electrons
                MM = qel.*[SeTJ1.*dRnTJ(in1) SeTJ2.*dRnTJ(in2) SeTJ3.*dRnTJ(in3)];
                Jmat0 = Jmat0 + sparse(ijs+nn,nn+IRBTJP(ijs),MM(mask_ijs),neq,neq);
                % Derivative of recombination term w.r.t. holes
                MM = qel.*[SeTJ1.*dRpTJ(in1) SeTJ2.*dRpTJ(in2) SeTJ3.*dRpTJ(in3)];
                Jmat0 = Jmat0 + sparse(ijs+nn,pp+ILBTJP(ijs),MM(mask_ijs),neq,neq);
                % Derivative of recombination term w.r.t. phi (left node)
                MM = qel.*[SeTJ1.*dRphin(in1) SeTJ2.*dRphin(in2) SeTJ3.*dRphin(in3)];
                Jmat0 = Jmat0 + sparse(ijs+nn,IRBTJP(ijs),MM(mask_ijs),neq,neq);
                % Derivative of recombination term w.r.t. phi (left node)
                MM = qel.*[SeTJ1.*dRphip(in1) SeTJ2.*dRphip(in2) SeTJ3.*dRphip(in3)];
                Jmat0 = Jmat0 + sparse(ijs+nn,ILBTJP(ijs),MM(mask_ijs),neq,neq);
                %
                % Assembling in hole continuity equation
                % Recombination term
                MM = qel.*[SeTJ1.*GBTBT(in1) SeTJ2.*GBTBT(in2) SeTJ3.*GBTBT(in3)];
                tvet = tvet + sparse(ijs+pp,1,MM(mask_ijs),neq,1);
                %
                % Derivative of recombination term w.r.t. holes
                MM = qel.*[SeTJ1.*dRpTJ(in1) SeTJ2.*dRpTJ(in2) SeTJ3.*dRpTJ(in3)];
                Jmat0 = Jmat0 + sparse(ijs+pp,pp+ILBTJP(ijs),MM(mask_ijs),neq,neq);
                % Derivative of recombination term w.r.t. electrons
                MM = qel.*[SeTJ1.*dRnTJ(in1) SeTJ2.*dRnTJ(in2) SeTJ3.*dRnTJ(in3)];
                Jmat0 = Jmat0 + sparse(ijs+pp,nn+IRBTJP(ijs),MM(mask_ijs),neq,neq);
                % Derivative of recombination term w.r.t. phi (left node)
                MM = qel.*[SeTJ1.*dRphin(in1) SeTJ2.*dRphin(in2) SeTJ3.*dRphin(in3)];
                Jmat0 = Jmat0 + sparse(ijs+pp,IRBTJP(ijs),MM(mask_ijs),neq,neq);
                % Derivative of recombination term w.r.t. phi (left node)
                MM = qel.*[SeTJ1.*dRphip(in1) SeTJ2.*dRphip(in2) SeTJ3.*dRphip(in3)];
                Jmat0 = Jmat0 + sparse(ijs+pp,ILBTJP(ijs),MM(mask_ijs),neq,neq);
            end
        end
    end
    
    %     MM = [Lp1.*Rrad2D(iiQW1) Lp2.*Rrad2D(iiQW2)];
    MM=[Se1.*Rrad(in1) Se2.*Rrad(in2) Se3.*Rrad(in3)];
    PspBulk=1000*sum(sparse(MM))*h*(Clight*1e-2/mean(1e-9*mode.vlambda)); % milliwatt
    
    mode.RSRH = RSRH;
    mode.Rrad = Rrad;
    mode.RAuger = RAuger;
    mode.GBTBT = GBTBT;
end
mode.PspBulk=PspBulk;

% =============================================================================================100
% assem Poisson equation
% =============================================================================================100
% #############################################################################################100
if((isfield(mode,'ionization'))&&(strcmp(mode.ionization,'Incomplete')))
    gD = 2; n1 = Nc.*exp(-mesh.DeltaEd./Vt);
    dop_dp  =    dop_d./( 1 + gD.*nF./(n1.*gamman)); % ionized donor density, 1/cm^3
    ddop_dp =  - dop_d./((1 + gD.*nF./(n1.*gamman)).^2) ...
        .*gD.*((n1.*gamman - nF.*n1.*dgamman)./((n1.*gamman).^2));
    dop_dp(dop_d==0)=0;
    ddop_dp(dop_d==0)=0;
    %
    gA = 4; p1 = Nv.*exp(-mesh.DeltaEa./Vt);
    dop_am  =    dop_a./( 1 + gA.*pF./(p1.*gammap)); % ionized donor density, 1/cm^3
    ddop_am =  - dop_a./((1 + gA.*pF./(p1.*gammap)).^2) ...
        .*gA.*((p1.*gammap - pF.*p1.*dgammap)./((p1.*gammap).^2));
    dop_am(dop_a==0)=0;
    ddop_am(dop_a==0)=0;
    %
else
    dop_dp = dop_d;
    dop_am = dop_a;
end
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
% =============================================================================================100
qelNorm=mode.qelNorm;
qelNorm2D=mode.qelNorm2D;
% doping charge in Poisson rhs
MM=qelNorm.*[Se1.*dop_dp(in1) Se2.*dop_dp(in2) Se3.*dop_dp(in3)];
tvet=tvet+sparse(ijr,1,-MM(mask_ijr),neq,1);
%
MM=qelNorm.*[Se1.*dop_am(in1) Se2.*dop_am(in2) Se3.*dop_am(in3)];
tvet=tvet+sparse(ijr,1,MM(mask_ijr),neq,1);
%
% incomplete ionization
if((isfield(mode,'ionization'))&&(strcmp(mode.ionization,'Incomplete')))
    if(mode.nflg)
        MM=qelNorm.*[ddop_dp(in1).*Se1 ddop_dp(in2).*Se2 ddop_dp(in3).*Se3];
        Jmat0=Jmat0+sparse(ijr,ijr+nn,-MM(mask_ijr),neq,neq);
    end
    if(mode.pflg)
        MM=qelNorm.*[ddop_am(in1).*Se1 ddop_am(in2).*Se2 ddop_am(in3).*Se3];
        Jmat0=Jmat0+sparse(ijr,ijr+pp,MM(mask_ijr),neq,neq);
    end
end
% =============================================================================================100
% elec charge in Poisson rhs
if(mode.nflg), MM=qelNorm.*[Se1 Se2 Se3];
    Kmat0 = Kmat0 + sparse(ijr,ijr+nn,MM(mask_ijr),neq,neq);
end
% =============================================================================================100
% hole charge in Poisson rhs
if(mode.pflg), MM=-qelNorm.*[Se1 Se2 Se3];
    Kmat0 = Kmat0 + sparse(ijr,ijr+pp,MM(mask_ijr),neq,neq);
end
% #############################################################################################100
% =============================================================================================100
% assem electron continuity equation
% =============================================================================================100
if(mode.nflg) % ##############################################################################100
    %
    
    
    VT12=(Vt(in1)+Vt(in2))/2;
    VT13=(Vt(in1)+Vt(in3))/2;
    VT23=(Vt(in2)+Vt(in3))/2;
    
    offset_n=zeros(1,nn);
    offset_nn=zeros(1,nn);
    
    % Bernoulli functions
    offset_n(iq)=affinity(iq);  % Vt(iq).*log(Nc(iq))
    offset_nn(iq)=log(Nc(iq)) ;
    offset_n1=offset_n(in1); offset_n2=offset_n(in2); offset_n3=offset_n(in3);
    offset_nn1=offset_nn(in1) ;offset_nn2=offset_nn(in2) ;offset_nn3=offset_nn(in3) ;
    
    
    
    delta12 = (phi1 + offset_n1 )./VT12+log(gamman(in1))+offset_nn1 - (phi2 + offset_n2 )./VT12-log(gamman(in2))-offset_nn2;
    delta23 = (phi2 + offset_n2 )./VT23+log(gamman(in2))+offset_nn2 - (phi3 + offset_n3 )./VT23-log(gamman(in3))-offset_nn3;
    delta31 = (phi3 + offset_n3 )./VT13+log(gamman(in3))+offset_nn3 - (phi1 + offset_n1 )./VT13-log(gamman(in1))-offset_nn1;
    %
    B12=G3.*bern(delta12);  dB12=G3.*dbern(delta12)./VT12;
    B21=G3.*bern(-delta12); dB21=G3.*dbern(-delta12)./VT12;
    B23=G1.*bern(delta23);  dB23=G1.*dbern(delta23)./VT23;
    B32=G1.*bern(-delta23); dB32=G1.*dbern(-delta23)./VT23;
    B31=G2.*bern(delta31);  dB31=G2.*dbern(delta31)./VT13;
    B13=G2.*bern(-delta31); dB13=G2.*dbern(-delta31)./VT13;
    % =============================================================================================100
    dDndphi1=dDndphi1.*Vt1; dDndphi2=dDndphi2.*Vt2; dDndphi3=dDndphi3.*Vt3;
    % =============================================================================================100
    M11=Dn.*(s3./l3.*B12+s2./l2.*B13);
    M12=-Dn.*s3./l3.*B21;
    M13=-Dn.*s2./l2.*B31;
    M21=-Dn.*s3./l3.*B12;
    M22=Dn.*(s1./l1.*B23+s3./l3.*B21);
    M23=-Dn.*s1./l1.*B32;
    M31=-Dn.*s2./l2.*B13;
    M32=-Dn.*s1./l1.*B23;
    M33=Dn.*(s2./l2.*B31+s1./l1.*B32);
    MM=qel.*[M11 M12 M13 M21 M22 M23 M31 M32 M33];
    %
    Kmat0=Kmat0+sparse(iis+nn,jjs+nn,MM(mask_iis),neq,neq); % Assembly electron eqs.
    Kmat0=Kmat0+sparse(iit   ,jjt+nn,-MM(mask_iit),neq,neq); % Assembly current eq.
    % derivatives w.r.t. potential ================================================================100
    M11=Dn.*((s3./l3.*dB12+s2./l2.*dB13).*elec1+s3./l3.*dB21.*elec2+s2./l2.*dB31.*elec3);
    M12=Dn.*(-s3./l3.*dB12.*elec1-s3./l3.*dB21.*elec2);
    M13=Dn.*(-s2./l2.*dB13.*elec1-s2./l2.*dB31.*elec3);
    M21=Dn.*(-s3./l3.*dB12.*elec1-s3./l3.*dB21.*elec2);
    M22=Dn.*(s3./l3.*dB12.*elec1+(s1./l1.*dB23+s3./l3.*dB21).*elec2+s1./l1.*dB32.*elec3);
    M23=Dn.*(-s1./l1.*dB23.*elec2-s1./l1.*dB32.*elec3);
    M31=Dn.*(-s2./l2.*dB13.*elec1-s2./l2.*dB31.*elec3);
    M32=Dn.*(-s1./l1.*dB23.*elec2-s1./l1.*dB32.*elec3);
    M33=Dn.*(s2./l2.*dB13.*elec1+s1./l1.*dB23.*elec2+(s2./l2.*dB31+s1./l1.*dB32).*elec3);
    MM=qel.*[M11 M12 M13 M21 M22 M23 M31 M32 M33];
    %
    Jmat0=Jmat0+sparse(iis+nn,jjs,MM(mask_iis),neq,neq); % Assembly electron eqs.
    Jmat0=Jmat0+sparse(iit   ,jjt,-MM(mask_iit),neq,neq); % Assembly current eq.
    % Fermi statistics ============================================================================100
    if((isfield(mode,'stats'))&&(strcmp(mode.stats,'Fermi')))
        
        M11=Dn.*((s3./l3.*VT12.*dB12+s2./l2.*VT13.*dB13).*elec1+s3./l3.*VT12.*dB21.*elec2+s2./l2.*VT13.*dB31.*elec3);
        M12=Dn.*(-s3./l3.*VT12.*dB12.*elec1-s3./l3.*VT12.*dB21.*elec2);
        M13=Dn.*(-s2./l2.*VT13.*dB13.*elec1-s2./l2.*VT13.*dB31.*elec3);
        M21=Dn.*(-s3./l3.*VT12.*dB12.*elec1-s3./l3.*VT12.*dB21.*elec2);
        M22=Dn.*(s3./l3.*VT12.*dB12.*elec1+(s1./l1.*VT23.*dB23+s3./l3.*VT12.*dB21).*elec2+s1./l1.*VT23.*dB32.*elec3);
        M23=Dn.*(-s1./l1.*VT23.*dB23.*elec2-s1./l1.*VT23.*dB32.*elec3);
        M31=Dn.*(-s2./l2.*VT13.*dB13.*elec1-s2./l2.*VT13.*dB31.*elec3);
        M32=Dn.*(-s1./l1.*VT23.*dB23.*elec2-s1./l1.*VT23.*dB32.*elec3);
        M33=Dn.*(s2./l2.*VT13.*dB13.*elec1+s1./l1.*VT23.*dB23.*elec2+(s2./l2.*VT13.*dB31+s1./l1.*VT23.*dB32).*elec3);
        MM=qel.*[M11.*dlGn1 M12.*dlGn2 M13.*dlGn3 ...
            M21.*dlGn1 M22.*dlGn2 M23.*dlGn3 ...
            M31.*dlGn1 M32.*dlGn2 M33.*dlGn3];
        
        %
        Jmat0=Jmat0+sparse(iis+nn,jjs+nn,MM(mask_iis),neq,neq); % Assembly electron eqs.
        Jmat0=Jmat0+sparse(iit   ,jjt+nn,-MM(mask_iit),neq,neq); % Assembly current eq.
    end
    
    % =============================================================================================100
    M11=((s3./l3.*B12+s2./l2.*B13).*elec1-s3./l3.*B21.*elec2-s2./l2.*B31.*elec3).*dDndphi1;
    M12=((s3./l3.*B12+s2./l2.*B13).*elec1-s3./l3.*B21.*elec2-s2./l2.*B31.*elec3).*dDndphi2;
    M13=((s3./l3.*B12+s2./l2.*B13).*elec1-s3./l3.*B21.*elec2-s2./l2.*B31.*elec3).*dDndphi3;
    M21=(-s3./l3.*B12.*elec1+(s1./l1.*B23+s3./l3.*B21).*elec2-s1./l1.*B32.*elec3).*dDndphi1;
    M22=(-s3./l3.*B12.*elec1+(s1./l1.*B23+s3./l3.*B21).*elec2-s1./l1.*B32.*elec3).*dDndphi2;
    M23=(-s3./l3.*B12.*elec1+(s1./l1.*B23+s3./l3.*B21).*elec2-s1./l1.*B32.*elec3).*dDndphi3;
    M31=(-s2./l2.*B13.*elec1-s1./l1.*B23.*elec2+(s2./l2.*B31+s1./l1.*B32).*elec3).*dDndphi1;
    M32=(-s2./l2.*B13.*elec1-s1./l1.*B23.*elec2+(s2./l2.*B31+s1./l1.*B32).*elec3).*dDndphi2;
    M33=(-s2./l2.*B13.*elec1-s1./l1.*B23.*elec2+(s2./l2.*B31+s1./l1.*B32).*elec3).*dDndphi3;
    MM=qel.*[M11 M12 M13 M21 M22 M23 M31 M32 M33];
    %
    Jmat0=Jmat0+sparse(iis+nn,jjs,MM(mask_iis),neq,neq); % Assembly electron eqs.
    Jmat0=Jmat0+sparse(iit   ,jjt,-MM(mask_iit),neq,neq); % Assembly current eq.
    % =============================================================================================100
    % direct recombination
    MM=qel.*[Se1.*R(in1) Se2.*R(in2) Se3.*R(in3)];
    tvet=tvet+sparse(ijs+nn,1,MM(mask_ijs),neq,1);
    % direct recombination, derivatives
    MM=qel.*[Se1.*dRn(in1) Se2.*dRn(in2) Se3.*dRn(in3)];
    Jmat0=Jmat0+sparse(ijs+nn,ijs+nn,MM(mask_ijs),neq,neq);
    if(mode.pflg)
        MM=qel.*[Se1.*dRp(in1) Se2.*dRp(in2) Se3.*dRp(in3)];
        Jmat0=Jmat0+sparse(ijs+nn,ijs+pp,MM(mask_ijs),neq,neq);
    end
    if(mode.AugerExpulsion3D==1)
        RAugerExpulsion=zeros(1,nn);
        dRAugerExpulsion_n3D=zeros(1,nn);
        dRAugerExpulsion_p3D=zeros(1,nn);
        RAugerExpulsion(mesh.ICAV)=RAuger_n3D(mesh.ICAV);
        dRAugerExpulsion_n3D(mesh.ICAV)=dRAuger_n3D_n3D(mesh.ICAV);
        dRAugerExpulsion_p3D(mesh.ICAV)=dRAuger_n3D_p3D(mesh.ICAV);
        %
        MM=qel.*[SeCAV1.*RAugerExpulsion(in1) SeCAV2.*RAugerExpulsion(in2) SeCAV3.*RAugerExpulsion(in3)];
        tvet=tvet+sparse(ijs+nn,1,MM(mask_ijs),neq,1);
        MM=qel.*[SeCAV1.*dRAugerExpulsion_n3D(in1) SeCAV2.*dRAugerExpulsion_n3D(in2) SeCAV3.*dRAugerExpulsion_n3D(in3)];
        Jmat0=Jmat0+sparse(ijs+nn,ijs+nn,MM(mask_ijs),neq,neq);
        if(mode.pflg)
            MM=qel.*[SeCAV1.*dRAugerExpulsion_p3D(in1) SeCAV2.*dRAugerExpulsion_p3D(in2) SeCAV3.*dRAugerExpulsion_p3D(in3)];
            Jmat0=Jmat0+sparse(ijs+nn,ijs+pp,MM(mask_ijs),neq,neq);
        end
    end
    % =============================================================================================100
    % Schottky contacts
    % =============================================================================================100
    for ic=1:nc
        ii=find((mesh.contact==ic) & (mesh.iq)); % semiconductor nodes on contact ic
        nq=length(ii);
        if((geom.contact_type(ic)==2) && (nq>0))
            jj=1:nq;
            ll(jj)=mesh.contact_seg(ic,ii);
            vsurfn(jj)=mesh.vsurfn(ii);
            ii=find((contact==ic) & iq);
            n_eq(1:nq)=nint(ii).*exp(phi(ii)-v_ls((ic-1)+1)./Vt(ii));
            tvet(ii+nn)=tvet(ii+nn)-(vsurfn.*n_eq.*ll).';
            Kmat0=Kmat0+sparse(ii+nn,ii+nn,vsurfn.*ll,neq,neq);
        end
    end
    % =============================================================================================100
end % #########################################################################################100
% =============================================================================================100
% assem hole continuity equation
% =============================================================================================100
if(mode.pflg) % ##############################################################################100
    %
    offset_p=zeros(1,nn);
    offset_pp=zeros(1,nn);
    
    offset_p(iq)=affinity(iq)+Eg(iq);
    offset_pp(iq)=-log(Nv(iq));
    offset_p1=offset_p(in1); offset_p2=offset_p(in2); offset_p3=offset_p(in3);
    offset_pp1=offset_pp(in1); offset_pp2=offset_pp(in2); offset_pp3=offset_pp(in3);
    
    % Bernoulli functions
    delta12 = (phi1 + offset_p1)./VT12-log(gammap(in1))+offset_pp1 - (phi2 + offset_p2)./VT12+log(gammap(in2))-offset_pp2;
    delta23 = (phi2 + offset_p2 )./VT23-log(gammap(in2))+offset_pp2 - (phi3 + offset_p3 )./VT23+log(gammap(in3))-offset_pp3;
    delta31 = (phi3 + offset_p3)./VT13-log(gammap(in3))+offset_pp3 - (phi1 + offset_p1 )./VT13+log(gammap(in1))-offset_pp1;
    %
    B12=G3.*bern(delta12);  dB12=G3.*dbern(delta12)./VT12;
    B21=G3.*bern(-delta12); dB21=G3.*dbern(-delta12)./VT12;
    B23=G1.*bern(delta23);  dB23=G1.*dbern(delta23)./VT23;
    B32=G1.*bern(-delta23); dB32=G1.*dbern(-delta23)./VT23;
    B31=G2.*bern(delta31);  dB31=G2.*dbern(delta31)./VT13;
    B13=G2.*bern(-delta31); dB13=G2.*dbern(-delta31)./VT13;
    % =============================================================================================100
    dDpdphi1=dDpdphi1.*Vt1; dDpdphi2=dDpdphi2.*Vt2; dDpdphi3=dDpdphi3.*Vt3;
    % =============================================================================================100
    M11=Dp.*(s3./l3.*B21+s2./l2.*B31);
    M12=-Dp.*s3./l3.*B12;
    M13=-Dp.*s2./l2.*B13;
    M21=-Dp.*s3./l3.*B21;
    M22=Dp.*(s1./l1.*B32+s3./l3.*B12);
    M23=-Dp.*s1./l1.*B23;
    M31=-Dp.*s2./l2.*B31;
    M32=-Dp.*s1./l1.*B32;
    M33=Dp.*(s2./l2.*B13+s1./l1.*B23);
    MM=qel.*[M11 M12 M13 M21 M22 M23 M31 M32 M33];
    %
    Kmat0=Kmat0+sparse(iis+pp,jjs+pp,MM(mask_iis),neq,neq); % Assembly hole eqs.
    Kmat0=Kmat0+sparse(iit   ,jjt+pp,MM(mask_iit),neq,neq); % Assembly current eq.
    % derivatives w.r.t. potential ================================================================100
    M11=Dp.*((-s3./l3.*dB21-s2./l2.*dB31).*hole1-s3./l3.*dB12.*hole2-s2./l2.*hole3.*dB13);
    M12=Dp.*(s3./l3.*dB21.*hole1+s3./l3.*dB12.*hole2);
    M13=Dp.*(s2./l2.*dB31.*hole1+s2./l2.*hole3.*dB13);
    M21=Dp.*(s3./l3.*dB21.*hole1+s3./l3.*dB12.*hole2);
    M22=Dp.*(-s3./l3.*dB21.*hole1+(-s1./l1.*dB32-s3./l3.*dB12).*hole2-s1./l1.*hole3.*dB23);
    M23=Dp.*(s1./l1.*dB32.*hole2+s1./l1.*hole3.*dB23);
    M31=Dp.*(s2./l2.*dB31.*hole1+s2./l2.*hole3.*dB13);
    M32=Dp.*(s1./l1.*dB32.*hole2+s1./l1.*hole3.*dB23);
    M33=Dp.*(-s2./l2.*dB31.*hole1-s1./l1.*dB32.*hole2+(-s2./l2.*dB13-s1./l1.*dB23).*hole3);
    MM=qel.*[M11 M12 M13 M21 M22 M23 M31 M32 M33];
    %
    Jmat0=Jmat0+sparse(iis+pp,jjs,MM(mask_iis),neq,neq);
    Jmat0=Jmat0+sparse(iit   ,jjt,MM(mask_iit),neq,neq);
    % Fermi statistics ============================================================================100
    if((isfield(mode,'stats'))&&(strcmp(mode.stats,'Fermi')))
        M11=Dp.*((-s3./l3.*dB21.*VT12-s2./l2.*dB31.*VT13).*hole1-s3./l3.*dB12.*VT12.*hole2-s2./l2.*hole3.*dB13.*VT13);
        M12=Dp.*(s3./l3.*dB21.*VT12.*hole1+s3./l3.*dB12.*VT12.*hole2);
        M13=Dp.*(s2./l2.*dB31.*VT13.*hole1+s2./l2.*hole3.*dB13.*VT13);
        M21=Dp.*(s3./l3.*dB21.*VT12.*hole1+s3./l3.*dB12.*VT12.*hole2);
        M22=Dp.*(-s3./l3.*dB21.*VT12.*hole1+(-s1./l1.*dB32.*VT23-s3./l3.*dB12.*VT12).*hole2-s1./l1.*hole3.*dB23.*VT23);
        M23=Dp.*(s1./l1.*dB32.*VT23.*hole2+s1./l1.*hole3.*dB23.*VT23);
        M31=Dp.*(s2./l2.*dB31.*VT13.*hole1+s2./l2.*hole3.*dB13.*VT13);
        M32=Dp.*(s1./l1.*dB32.*VT23.*hole2+s1./l1.*hole3.*dB23.*VT23);
        M33=Dp.*(-s2./l2.*dB31.*VT13.*hole1-s1./l1.*dB32.*VT23.*hole2+(-s2./l2.*dB13.*VT13-s1./l1.*dB23.*VT23).*hole3);
        MM=-qel.*[M11.*dlGp1 M12.*dlGp2 M13.*dlGp3 ...
            M21.*dlGp1 M22.*dlGp2 M23.*dlGp3 ...
            M31.*dlGp1 M32.*dlGp2 M33.*dlGp3];
        %
        Jmat0=Jmat0+sparse(iis+pp,jjs+pp,MM(mask_iis),neq,neq); % Assembly hole eqs.
        Jmat0=Jmat0+sparse(iit   ,jjt+pp,MM(mask_iit),neq,neq); % Assembly current eq.
    end
    % =============================================================================================100
    M11=((s3./l3.*B21+s2./l2.*B31).*hole1-s3./l3.*B12.*hole2-s2./l2.*hole3.*B13).*dDpdphi1;
    M12=((s3./l3.*B21+s2./l2.*B31).*hole1-s3./l3.*B12.*hole2-s2./l2.*hole3.*B13).*dDpdphi2;
    M13=((s3./l3.*B21+s2./l2.*B31).*hole1-s3./l3.*B12.*hole2-s2./l2.*hole3.*B13).*dDpdphi3;
    M21=(-s3./l3.*B21.*hole1+(s1./l1.*B32+s3./l3.*B12).*hole2-s1./l1.*hole3.*B23).*dDpdphi1;
    M22=(-s3./l3.*B21.*hole1+(s1./l1.*B32+s3./l3.*B12).*hole2-s1./l1.*hole3.*B23).*dDpdphi2;
    M23=(-s3./l3.*B21.*hole1+(s1./l1.*B32+s3./l3.*B12).*hole2-s1./l1.*hole3.*B23).*dDpdphi3;
    M31=(-s2./l2.*B31.*hole1-s1./l1.*B32.*hole2+(s2./l2.*B13+s1./l1.*B23).*hole3).*dDpdphi1;
    M32=(-s2./l2.*B31.*hole1-s1./l1.*B32.*hole2+(s2./l2.*B13+s1./l1.*B23).*hole3).*dDpdphi2;
    M33=(-s2./l2.*B31.*hole1-s1./l1.*B32.*hole2+(s2./l2.*B13+s1./l1.*B23).*hole3).*dDpdphi3;
    MM=qel.*[M11 M12 M13 M21 M22 M23 M31 M32 M33];
    %
    Jmat0=Jmat0+sparse(iis+pp,jjs,MM(mask_iis),neq,neq); % Assembly hole eqs.
    Jmat0=Jmat0+sparse(iit   ,jjt,MM(mask_iit),neq,neq); % Assembly current eq.
    % =============================================================================================100
    % direct recombination, derivatives
    MM=qel.*[Se1.*R(in1) Se2.*R(in2) Se3.*R(in3)];
    tvet=tvet+sparse(ijs+pp,1,MM(mask_ijs),neq,1);
    % direct recombination, derivatives
    MM=qel.*[Se1.*dRp(in1) Se2.*dRp(in2) Se3.*dRp(in3)];
    Jmat0=Jmat0+sparse(ijs+pp,ijs+pp,MM(mask_ijs),neq,neq);
    if(mode.nflg)
        MM=qel.*[Se1.*dRn(in1) Se2.*dRn(in2) Se3.*dRn(in3)];
        Jmat0=Jmat0+sparse(ijs+pp,ijs+nn,MM(mask_ijs),neq,neq);
    end
    if(mode.AugerExpulsion3D==1)
        RAugerExpulsion=zeros(1,nn);
        dRAugerExpulsion_p3D=zeros(1,nn);
        dRAugerExpulsion_n3D=zeros(1,nn);
        RAugerExpulsion(mesh.ICAV)=RAuger_p3D(mesh.ICAV);
        dRAugerExpulsion_p3D(mesh.ICAV)=dRAuger_p3D_p3D(mesh.ICAV);
        dRAugerExpulsion_n3D(mesh.ICAV)=dRAuger_p3D_n3D(mesh.ICAV);
        %
        MM=qel.*[SeCAV1.*RAugerExpulsion(in1) SeCAV2.*RAugerExpulsion(in2) SeCAV3.*RAugerExpulsion(in3)];
        tvet=tvet+sparse(ijs+pp,1,MM(mask_ijs),neq,1);
        MM=qel.*[SeCAV1.*dRAugerExpulsion_p3D(in1) SeCAV2.*dRAugerExpulsion_p3D(in2) SeCAV3.*dRAugerExpulsion_p3D(in3)];
        Jmat0=Jmat0+sparse(ijs+pp,ijs+pp,MM(mask_ijs),neq,neq);
        if(mode.nflg)
            MM=qel.*[SeCAV1.*dRAugerExpulsion_n3D(in1) SeCAV2.*dRAugerExpulsion_n3D(in2) SeCAV3.*dRAugerExpulsion_n3D(in3)];
            Jmat0=Jmat0+sparse(ijs+pp,ijs+nn,MM(mask_ijs),neq,neq);
        end
    end
    % =============================================================================================100
    % Schottky contacts
    for ic=1:nc
        ii=find((mesh.contact==ic)&(mesh.iq)); % semiconductor nodes on contact ic
        nq=length(ii);
        if((geom.contact_type(ic)==2)&&(nq>0))
            ll(1:nq) = mesh.contact_seg(ic,ii);
            vsurfp(1:nq) = mesh.vsurfp(ii);
            ii=find((contact==ic)&iq);
            jj=ii;
            p_eq(1:nq) = nint(ii).*exp(-phi(ii)+v_ls((ic-1)+1)./Vt(ii));
            tvet(ii+pp)=tvet(ii+pp)-(vsurfp.*p_eq.*ll).';
            Kmat0=Kmat0+sparse(ii+pp,ii+pp,vsurfp.*ll,neq,neq);
        end
    end
    % =============================================================================================100
end % #########################################################################################100
% =============================================================================================100
% =============================================================================================100
% Add circuit equations: V-Z*I-V0=0 (V-driven); I-Y*V-I0=0 (I-driven)
% =============================================================================================100
Zmat=sparse(mode.Zmat);
Jmat2=sparse(neq,neq);

if mode.IdriveON==1   % current driving
    Jmat2=Jmat2+sparse(qq+1,rr+1,1,neq,neq); % I
    %
    for ic=1:nm
        ii=(ic-1)+1;
        for jc=1:nm
            jj=(jc-1)+1;
            Jmat2(qq+ii,qq+jj)=Ymat(ii,jj);
        end
    end % -Y*V
    tvet(qq+(1:nm))=-i0_dd; % -V0
else    % voltage driving
    Jmat2=Jmat2+sparse(qq+1,qq+1,1,neq,neq); % V
    %
    for ic=1:nm
        ii=(ic-1)+1;
        for jc=1:nm
            jj=(jc-1)+1;
            Jmat2(qq+ii,rr+jj)=Zmat(ii,jj);
        end
    end % -Z*I
    tvet(qq+(1:nm))=-v0_dd; % -V0
end
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
if(mode.oflg)
    iTappo=mode.iTappo;
    N2Dtot=0; P2Dtot=0; Psp=0;
    mode.gmod=zeros(nmodes,1);
    mode.matgain=0;
    %        mode.DeltaN=0;
    mode.RSRHQW=zeros(1,mesh.nn);
    mode.RradQW=zeros(1,mesh.nn);
    mode.RAugerQW=zeros(1,mesh.nn);
    mode.RLeakageQW=zeros(1,mesh.nn);
    mode.RradQW_SBE=zeros(1,mesh.nn);
    mode.Ccapn3D=zeros(1,mesh.nn);
    mode.Ccapp3D=zeros(1,mesh.nn);
    IntCcapp=[0 0];
    IntCcapn=[0 0];
    GaSum=0;
    IntSRH=0;
    IntAug=0;
    IntStim=0;
    IntRad=0;
    IntLea=0;
    IntRecQW3Dn=0; % integral over all the QW (3D) recombinations
    IntRecQW3Dp=0; % integral over all the QW (3D) recombinations
    ImatQW = sparse(neq,neq);
    tvetQW = sparse(neq,1);
    MMgain=0;
    for indQW=1:NQW
        %
        % Basic parameters: geometry
        WQW = mesh.vWMQW{indQW}; %-- cm
        IQW = unique(mesh.IMQW{indQW}); % questo unique metterlo in rectmesh
        inQW = mesh.inMQW{indQW}; % index to QW nodes
        inQWP = mesh.inMQWP{indQW}; % including passivation
        VV = vv+(indQW-1)*length(mesh.xgrid);
        WW = ww+(indQW-1)*length(mesh.xgrid);
        %
        iiQW1=1:nnQW-1; iiQW2=2:nnQW;
        maskQW=true(1,nnQW);
        % maskQW([1,end])=0;
        % maskQW([1,2,end-1,end])=0;
        %
        SeQW1=zeros(1,nt); SeQW2=SeQW1; SeQW3=SeQW1;
        SeQW1(mesh.ITrMQW{indQW})=Se1(mesh.ITrMQW{indQW});
        SeQW2(mesh.ITrMQW{indQW})=Se2(mesh.ITrMQW{indQW});
        SeQW3(mesh.ITrMQW{indQW})=Se3(mesh.ITrMQW{indQW});
        ijrQW=[iiQW1 iiQW2]; inQW1=inQW(iiQW1); inQW2=inQW(iiQW2);
        %
        IIQWP = repmat(1:mesh.nnx,mesh.nny,1); IIQWP=IIQWP(:).'; % for Poisson
        
        INQWP=inQWP(IIQWP); % for Ccap recombinations in 3D continuity equations
        %
        xQW = node(1,inQW); yQW = node(2,inQW);
        xcQW = (xQW(2:end) + xQW(1:end-1))/2; % edge centers, cm
        ycQW = (yQW(2:end) + yQW(1:end-1))/2; % edge centers, cm
        LeQW = xQW(2:end) - xQW(1:end-1); % edge length, cm
        Lp = zeros(1,nnQW); % box length, cm
        Lp(1:(nnQW-1)) = LeQW/2; Lp(2:nnQW) = Lp(2:nnQW) + LeQW/2;
        %
        Lp1=Lp(iiQW1)/2; Lp1(1)=2*Lp1(1);
        Lp2=Lp(iiQW2)/2; Lp2(end)=2*Lp2(end);
        Gcyl=1;
        if((isfield(mode,'symmetry'))&&(strcmp(mode.symmetry,'Cylindrical-Y')))
            Lp1=2.*pi.*Lp1.*xcQW;
            Lp2=2.*pi.*Lp2.*xcQW;
            Gcyl=2.*pi.*xcQW;
            yQW=yQW(1);
        end
        %
        if((isfield(mode,'symmetry'))&&(strcmp(mode.symmetry,'Cylindrical-X')))
            Lp1=2.*pi.*Lp1.*ycQW;
            Lp2=2.*pi.*Lp2.*ycQW;
            Gcyl=2.*pi.*ycQW;
            xQW=xQW(1);
        end
        % % debug: check Lp1,Lp2: integral of the constant function on the domain
        % %        this is expected to be equal to the length (z) of the QW
        % sum(diag(sparse([inQW1 inQW2],[inQW1 inQW2],[Lp1 Lp2],nn,nn)),1)
        %
        % Basic parameters: physics
        s_LoadConstants
        Vt2D = Vt(inQW);
        
        if mode.TQWFake==1 && mode.indv> 2
            load('DTfit') ;
            Iq=mode.ii_dd(end)*1e3 ;
            TQWfit=interp1(Ifit,DTfit,Iq)+293 ;
            Vt2D = ones(1,length(inQW))*TQWfit*kB/qel ;
            disp('Warning! Tqw fake!')
        end
        
        %
        tauscatn = mesh.tauscatnMQW{indQW}/mode.FatAdiab; % capture time, s
        tauscatp = mesh.tauscatpMQW{indQW}/mode.FatAdiab; % capture time, s
        meffn2D = mesh.meffnMQW{indQW};
        meffp2D = mesh.meffpMQW{indQW};
        DeltaEc = mesh.DeltaEcQW{indQW};
        DeltaEv = mesh.DeltaEvQW{indQW};
        lambda_C=mesh.lambda_C{indQW};
        lambda_V=mesh.lambda_V{indQW};
        Psi_C2=mesh.Psi_C2{indQW};
        Psi_V2=mesh.Psi_V2{indQW};
        NBound_C=length(lambda_C); % number of conduction band bound states
        NBound_V=length(lambda_V); % number of valence band bound states
        Ec0 = ones(NBound_C,1)*mode.ecb(inQW);
        Ev0 = ones(NBound_V,1)*mode.evb(inQW);
        Em = Ec0 - DeltaEc + lambda_C*ones(1,nnQW); % eigenvalue referred to bottom of QW
        Ed = Ev0 + DeltaEv - lambda_V*ones(1,nnQW); % eigenvalue referred to top of QB
        %        'TAROCCO N2 @@@@@@@@'
        %        N2 =.85* 4*pi/h^2*m0*qel*(sum((meffn2D*ones(1,mesh.nnxQW)).*(Ec0-Em),1)/10000);
        if isfield(mode,'Fat_N2')
            Fat_N2=mode.Fat_N2;
            Fat_P2=mode.Fat_P2;
        else
            Fat_N2=1;
            Fat_P2=1;
        end
        
        N2 = Fat_N2 * 4*pi/h^2*m0*qel*(sum((meffn2D*ones(1,mesh.nnxQW{indQW})).*(Ec0-Em),1)/10000)/mode.CarrierNorm2D;
        P2 = Fat_P2 * 4*pi/h^2*m0*qel*(sum((meffp2D*ones(1,mesh.nnxQW{indQW})).*(Ed-Ev0),1)/10000)/mode.CarrierNorm2D;
        
        Nc2D = mesh.Nc2D{indQW};%*ones(1,nnQW);
        Nv2D = mesh.Nv2D{indQW};%*ones(1,nnQW);
        MVt2D_C=ones(NBound_C,1)*Vt2D;
        MVt2D_V=ones(NBound_V,1)*Vt2D;
        %
        %           Dn2Dd = mode.FAT_idiffusioneQW.*mesh.mobnQW{indQW}.*Vt2D(iiQW1); % 2D diffusivity
        %           Dp2Dd = mode.FAT_idiffusioneQW.*mesh.mobpQW{indQW}.*Vt2D(iiQW1);  % 2D diffusivity
        Dn2Dd = mode.FAT_idiffusioneQW_E.*mesh.mobnQW{indQW}.*Vt2D(iiQW1); % 2D diffusivity
        Dp2Dd = mode.FAT_idiffusioneQW_H.*mesh.mobpQW{indQW}.*Vt2D(iiQW1);  % 2D diffusivity
        
        Damb=1./(1./Dn2Dd+1./Dp2Dd);
        %
        if(mode.iambipolarQW==1)
            Dn2D = Damb;
            Dp2D = Damb;
        else
            Dn2D = Dn2Dd;
            Dp2D = Dp2Dd;
        end
        %
        EFn3D = ones(NBound_C,1)*mode.EFn(inQW);
        EFp3D = ones(NBound_V,1)*mode.EFp(inQW);
        %
        iivQW=[iiQW1 iiQW1 iiQW2 iiQW2]; jjvQW=[iiQW1 iiQW2 iiQW1 iiQW2];
        invQW=[inQW1 inQW1 inQW2 inQW2]; jnvQW=[inQW1 inQW2 inQW1 inQW2];
        iivQWcur=[iiQW1 iiQW1]; jjvQWcur=[iiQW1 iiQW2];
        mask_iivQW=maskQW(iivQW);
        mask_iivQWcur=maskQW(iivQWcur);
        %
        n2D=n2Dc{indQW}; p2D=p2Dc{indQW};
        n2Di=mode.n2Di{indQW}; p2Di=mode.p2Di{indQW};
        %
        % checking if n2D>N2 or p2D>P2
        % indNaN=find(isnan(n2D));
        % if(not(isempty(indNaN))),n2D(indNaN)=mean(n2D); end
        if(iTappo==1 | iTappo==3 | iTappo==6)
            indN2=find(n2D>=N2);if(not(isempty(indN2))),n2D(indN2)=0.99*N2(indN2);end
            % indNaN=find(isnan(p2D));
            % if(not(isempty(indNaN))),p2D(indNaN)=mean(p2D); end
            indP2=find(p2D>=P2);if(not(isempty(indP2))),p2D(indP2)=0.99*P2(indP2);end
        end
        %
        % Evaluation of 2D quasi-Fermi levels
        % 2D electrons
        
        if(iTappo==1 | iTappo==3 | iTappo==6 | iTappo==7)
            EFn2D=invferdr_n2D(n2D,Em,Ec0,Nc2D,MVt2D_C,EFp3D(1,:)-3,EFn3D(1,:)+3);
        elseif(iTappo==0 | iTappo==2 | iTappo==4  | iTappo==8)
            EFn2D=invferdr_n2D_0(n2D,Em,Ec0,Nc2D,MVt2D_C,EFp3D(1,:)-3,EFn3D(1,:)+3);
        end
        
        xnum=(ones(NBound_C,1)*EFn2D-Em)./MVt2D_C; xden=(ones(NBound_C,1)*EFn2D-Ec0)./MVt2D_C;
        if(iTappo==1 | iTappo==3 | iTappo==6 | iTappo==7)
            [n2Dm,dn2D_2] = ferdr2D(xnum,xden);
        elseif(iTappo==0 | iTappo==2 | iTappo==4  | iTappo==8)
            [n2Dm,dn2D_2] = ferdr2D_0(xnum,xden);
        end
        n2Dm=Nc2D.*n2Dm;
        dn2Dm_EFn2D=Nc2D.*dn2D_2./MVt2D_C;
        dn2D_EFn2D=sum(dn2Dm_EFn2D,1);
        %
        % 2D holes
        if(iTappo==1 | iTappo==3  | iTappo==6 | iTappo==7)
            EFp2D=invferdr_p2D(p2D,Ed,Ev0,Nv2D,MVt2D_V,EFp3D(1,:)-3,EFn3D(1,:)+3);
        elseif(iTappo==0 | iTappo==2 | iTappo==4  | iTappo==8)
            EFp2D=invferdr_p2D_0(p2D,Ed,Ev0,Nv2D,MVt2D_V,EFp3D(1,:)-3,EFn3D(1,:)+3);
        end
        xnum=(Ed-ones(NBound_V,1)*EFp2D)./MVt2D_V; xden=(Ev0-ones(NBound_V,1)*EFp2D)./MVt2D_V;
        if(iTappo==1 | iTappo==3  | iTappo==6 | iTappo==7)
            [p2Dm,dp2D_2] = ferdr2D(xnum,xden);
        elseif(iTappo==0 | iTappo==2 | iTappo==4  | iTappo==8)
            [p2Dm,dp2D_2] = ferdr2D_0(xnum,xden);
        end
        p2Dm=Nv2D.*p2Dm;
        dp2Dm_EFp2D=-Nv2D.*dp2D_2./MVt2D_V;
        dp2D_EFp2D=sum(dp2Dm_EFp2D,1);
        %
        % Assembing 2D electrons charge in Poisson equation
        n2DmPois = zeros(size(Nc2D,1),length(mesh.xgrid));
        dn2Dm_EFn2DPois = n2DmPois;
        n2DmPois(:,1:size(n2Dm,2))=n2Dm;
        dn2Dm_EFn2DPois(:,1:size(dn2Dm_EFn2D,2))=dn2Dm_EFn2D;
        N2D = reshape(Psi_C2.'*n2DmPois,1,nn);
        N2Dtot = N2Dtot + N2D;
        dN2D_EFn2D = reshape(Psi_C2.'*dn2Dm_EFn2DPois,1,nn);
        dn2D_EFn2Dvec = reshape(ones(NBound_C,mesh.nny).'*dn2Dm_EFn2DPois,1,nn);
        dN2D_n2D = dN2D_EFn2D./dn2D_EFn2Dvec;
        dN2D_n2D(isnan(dN2D_n2D)) = 0;
        % qelNorm (not qelNorm2D) because N2D is multiplied by Psi
        % (1/cm)
        MM = qelNorm2D.*[Se1.*N2D(in1) Se2.*N2D(in2) Se3.*N2D(in3)];
        tvet = tvet + sparse(ijr,1,MM(mask_ijr),neq,1);
        MM = qelNorm2D.*[Se1.*dN2D_n2D(in1) Se2.*dN2D_n2D(in2) Se3.*dN2D_n2D(in3)];
        Jmat0 = Jmat0 + sparse(ijr,VV+IIQWP(ijr),MM(mask_ijr),neq,neq);
        %
        % Assembing 2D holes charge in Poisson equation
        p2DmPois = zeros(size(Nv2D,1),length(mesh.xgrid));
        dp2Dm_EFp2DPois = p2DmPois;
        p2DmPois(:,1:size(p2Dm,2))=p2Dm;
        dp2Dm_EFp2DPois(:,1:size(dp2Dm_EFp2D,2))=dp2Dm_EFp2D;
        P2D = reshape(Psi_V2.'*p2DmPois,1,nn);
        P2Dtot = P2Dtot + P2D;
        dP2D_EFp2D = reshape(Psi_V2.'*dp2Dm_EFp2DPois,1,nn);
        dp2D_EFp2Dvec = reshape(ones(NBound_V,mesh.nny).'*dp2Dm_EFp2DPois,1,nn);
        dP2D_p2D = dP2D_EFp2D./dp2D_EFp2Dvec;
        dP2D_p2D(isnan(dP2D_p2D)) = 0;
        MM = - qelNorm2D.*[Se1.*P2D(in1) Se2.*P2D(in2) Se3.*P2D(in3)];
        tvet = tvet + sparse(ijr,1,MM(mask_ijr),neq,1);
        MM = - qelNorm2D.*[Se1.*dP2D_p2D(in1) Se2.*dP2D_p2D(in2) Se3.*dP2D_p2D(in3)];
        Jmat0 = Jmat0 + sparse(ijr,WW+IIQWP(ijr),MM(mask_ijr),neq,neq);
        %
        % 2-D electron continuity equation ----------------------------------------
        if mode.idiffusioneQW==2 | mode.idiffusioneQW==3
            if(mode.idiffusioneQW==3)
                EwC=Ec0(1,:)-DeltaEc;
                nB2D=mean(Nc2D,1).*exp((EFn2D-EwC)./Vt2D);
                % nB2D=sum(Nc2D,1).*exp((EFn2D-EwC)./Vt2D);
                % nB2D=Nc2D(1).*exp((EFn2D-EwC)./Vt2D);
                gamman2D=n2D./nB2D;
                lnGn2D=log(gamman2D);
                dlnGn2D_n2D=((1-n2D./(Vt2D.*dn2D_EFn2D))./nB2D)./gamman2D;
                delta12 = (phi(inQW1)+Vt2D(iiQW1).*lnGn2D(iiQW1))./Vt2D(iiQW1) - (phi(inQW2)+Vt2D(iiQW2).*lnGn2D(iiQW2))./Vt2D(iiQW2);
            elseif(mode.idiffusioneQW==2)
                delta12 = phi(inQW1)./Vt2D(iiQW1) - phi(inQW2)./Vt2D(iiQW2);
            end
            n2D1=n2D(iiQW1); n2D2=n2D(iiQW2);
            B12=bern(delta12); dB12=dbern(delta12)./Vt2D(iiQW1);
            B21=bern(-delta12); dB21=dbern(-delta12)./Vt2D(iiQW2);
        elseif mode.idiffusioneQW==1
            B12 = 1; B21 = 1; % diffusion-only equations
        end
        %
        if mode.idiffusioneQW~=0
            M11 =   Gcyl.*Dn2D./LeQW.*B12;
            M12 = - Gcyl.*Dn2D./LeQW.*B21;
            M21 = - Gcyl.*Dn2D./LeQW.*B12;
            M22 =   Gcyl.*Dn2D./LeQW.*B21;
            MM  =   qel.*[M11 M12 M21 M22];
            Kmat0 = Kmat0 + sparse(VV+iivQW(mask_iivQW),VV+jjvQW(mask_iivQW),MM(mask_iivQW),neq,neq);
            MM  =   qel.*[M11 M12];
            ImatQW = ImatQW + sparse(VV+iivQWcur(mask_iivQWcur),VV+jjvQWcur(mask_iivQWcur),-MM(mask_iivQWcur),neq,neq);
        end
        %
        if mode.idiffusioneQW==2 | mode.idiffusioneQW==3
            % Derivative w.r.t. potential
            M11 =   Gcyl.*Dn2D./LeQW.*(dB12.*n2D1+dB21.*n2D2);
            M12 = - Gcyl.*Dn2D./LeQW.*(dB21.*n2D2+dB12.*n2D1);
            M21 = - Gcyl.*Dn2D./LeQW.*(dB21.*n2D2+dB12.*n2D1);
            M22 =   Gcyl.*Dn2D./LeQW.*(dB21.*n2D2+dB12.*n2D1);
            MM  =   qel.*[M11 M12 M21 M22];
            Jmat0=Jmat0+sparse(VV+iivQW(mask_iivQW),jnvQW(mask_iivQW),MM(mask_iivQW),neq,neq);
            %
            % Derivative Fermi factor
            if(mode.idiffusioneQW==3)
                MM  =   qel.*[Vt2D(iiQW1).*M11.*dlnGn2D_n2D(iiQW1) ...
                    Vt2D(iiQW2).*M12.*dlnGn2D_n2D(iiQW2) ...
                    Vt2D(iiQW1).*M21.*dlnGn2D_n2D(iiQW1) ...
                    Vt2D(iiQW2).*M22.*dlnGn2D_n2D(iiQW2)];
                Jmat0=Jmat0+sparse(VV+iivQW(mask_iivQW),VV+jjvQW(mask_iivQW),MM(mask_iivQW),neq,neq);
            end
        end
        
        %
        % 2-D hole continuity equation --------------------------------------------
        if mode.idiffusioneQW==2 | mode.idiffusioneQW==3
            if(mode.idiffusioneQW==3)
                EwV=Ev0(1,:)+DeltaEv;
                pB2D=mean(Nv2D,1).*exp((EwV-EFp2D)./Vt2D); % arbitrary :-(
                % pB2D=sum(Nv2D,1).*exp((EwV-EFp2D)./Vt2D); % arbitrary :-(
                % pB2D=Nv2D(1).*exp((EwV-EFp2D)./Vt2D); % arbitrary :-(
                gammap2D=p2D./pB2D;
                lnGp2D=log(gammap2D);
                dlnGp2D_p2D=((1+p2D./(Vt2D.*dp2D_EFp2D))./pB2D)./gammap2D;
                delta12 = (phi(inQW1)-Vt2D(iiQW1).*lnGp2D(iiQW1))./Vt2D(iiQW1) - (phi(inQW2)-Vt2D(iiQW2).*lnGp2D(iiQW2))./Vt2D(iiQW2);
                %
            elseif(mode.idiffusioneQW==2)
                
                delta12 = phi(inQW1)./Vt2D(iiQW1) - phi(inQW2)./Vt2D(iiQW2);
            end
            p2D1=p2D(iiQW1); p2D2=p2D(iiQW2);
            B12=bern(delta12); dB12=dbern(delta12)./Vt2D(iiQW2);
            B21=bern(-delta12); dB21=dbern(-delta12)./Vt2D(iiQW1);
        elseif mode.idiffusioneQW==1
            B12 = 1; B21 = 1;  % diffusion-only equations
        end
        %
        if mode.idiffusioneQW~=0
            M11 =   Gcyl.*Dp2D./LeQW.*B21;
            M12 = - Gcyl.*Dp2D./LeQW.*B12;
            M21 = - Gcyl.*Dp2D./LeQW.*B21;
            M22 =   Gcyl.*Dp2D./LeQW.*B12;
            MM  =   qel.*[M11 M12 M21 M22];
            Kmat0 = Kmat0 + sparse(WW+iivQW(mask_iivQW),WW+jjvQW(mask_iivQW),MM(mask_iivQW),neq,neq);
            MM  =   qel.*[M11 M12];
            ImatQW = ImatQW + sparse(WW+iivQWcur(mask_iivQWcur),WW+jjvQWcur(mask_iivQWcur),MM(mask_iivQWcur),neq,neq);
        end
        %
        % Derivative w.r.t. potential
        if mode.idiffusioneQW==2 | mode.idiffusioneQW==3
            M11 =   Gcyl.*Dp2D./LeQW.*(-dB21.*p2D1-dB12.*p2D2);
            M12 = - Gcyl.*Dp2D./LeQW.*(-dB12.*p2D2-dB21.*p2D1);
            M21 = - Gcyl.*Dp2D./LeQW.*(-dB21.*p2D1-dB12.*p2D2);
            M22 =   Gcyl.*Dp2D./LeQW.*(-dB12.*p2D2-dB21.*p2D1);
            MM  =   qel.*[M11 M12 M21 M22];
            Jmat0 = Jmat0+sparse(WW+iivQW(mask_iivQW),jnvQW(mask_iivQW),MM(mask_iivQW),neq,neq);
            %
            % Derivative Fermi factor
            if(mode.idiffusioneQW==3)
                MM  = - qel.*[Vt2D(iiQW1).*M11.*dlnGp2D_p2D(iiQW1) ...
                    Vt2D(iiQW2).*M12.*dlnGp2D_p2D(iiQW2)
                    Vt2D(iiQW1).*M21.*dlnGp2D_p2D(iiQW1) ...
                    Vt2D(iiQW2).*M22.*dlnGp2D_p2D(iiQW2)];
                Jmat0 = Jmat0+sparse(WW+iivQW(mask_iivQW),WW+jjvQW(mask_iivQW),MM(mask_iivQW),neq,neq);
            end
        end
        %
        % 2D capture terms
        % 2D electrons capture terms
        %        tauscatnT=tauscatn.*exp(-(TQW-mode.T0)/mode.T_tauscat);
        %        tauscatpT=tauscatp.*exp(-(TQW-mode.T0)/mode.T_tauscat);
        tauscatnT=tauscatn.*ones(1,nnQW);
        tauscatpT=tauscatp.*ones(1,nnQW);
        
        if(mode.taucarrier~=1) % no Vallone
            if isfield(mode,'T_tauscat')
                TQW=mesh.T(inQW);
                T_tauscat=mode.T_tauscat;
                Tcap_EXP=mode.Tcap_EXP;
                
                if T_tauscat~=0
                    %                    TQWm=mean(TQW);
                    %                    tauscatnT=tauscatn.*exp((TQW-mode.T0)/mode.T_tauscat);
                    %                    tauscatpT=tauscatp.*exp((TQW-mode.T0)/mode.T_tauscat);
                    
                    
                    %                        if Tcap_EXP==Inf
                    %
                    %                            Fat_cap=exp((TQW-mode.T0)/mode.T_tauscat);
                    %                        elseif Tcap_EXP==-Inf
                    %
                    %                            DeT=mode.T_tauscat;
                    %                            tauRat=mode.tauRat;
                    %
                    %                        else
                    %                     Fat_cap=(1+(TQW-mode.T0)/mode.T_tauscat).^Tcap_EXP;
                    Fat_cap=exp((TQW-mode.T300)/mode.T_tauscat).^Tcap_EXP;
                    %                        end
                    
                    if mode.leakage>0
                        Tstart=mode.Tstart;
                        taubottom=mode.tL_high;
                        DeT=mode.leakage;
                        tautop=mode.tL_low;
                        Fat_capleak=f_RaisedCosine(TQW,Tstart,DeT,taubottom,tautop);
                    end
                    Fcap=1;
                    %                    tauscatnT=tauscatn.*Fcap;
                    %                    tauscatpT=tauscatp.*Fcap;
                    tauscatnT=tauscatn.*Fat_cap;
                    tauscatpT=tauscatp.*Fat_cap;
                    %
                    % Se decidiamo di usare il modello di Vallone
                    
                    %                    tauscatnT=tauscatn.*ones(size(Fat_cap));
                    %                    tauscatpT=tauscatp.*ones(size(Fat_cap));
                end
                tauRat=mode.tauRat;
                tauescn=tauscatn./Fat_cap*tauRat.*ones(1,nnQW);
                tauescp=tauscatp./Fat_cap*tauRat.*ones(1,nnQW);
                mode.FatCap=Fat_cap;
                
            end
            
        else % Vallone
            if(isfield(mode,'Fat_cap_e'))
                tauscatnT=tauscatn.*mode.Fat_cap_e;
                tauscatpT=tauscatp.*mode.Fat_cap_h;
            end
            tauRat=mode.tauRat;
            tauescn=tauscatn*tauRat.*ones(1,nnQW);
            tauescp=tauscatp*tauRat.*ones(1,nnQW);
            mode.FatCap=mode.Fat_cap_e;
            
        end
        
        tauleakn=tauescn;
        tauleakp=tauescp;
        %                    tauescn=tauscatn.*tauRat.*ones(size(Fat_cap));
        %                    tauescp=tauscatp.*tauRat.*ones(size(Fat_cap));
        %                    tauescp=tauscatn.*Fat_cap
        
        dN3D_Efn3D =    dnF(inQW)./Vt2D;
        
        if(iTappo==1) % esattamente il termine degli svizzeri
            FCapn=exp((EFn2D-EFn3D(1,:))./Vt2D);
            Ccapn = - (1-FCapn).*(1-n2D./N2)./tauscatnT.*nF(inQW);
            dCcapn_N3D =  - (1-FCapn).*(1-n2D./N2)./tauscatnT + ...
                - FCapn.*(1-n2D./N2).*nF(inQW)./tauscatnT./dN3D_Efn3D./Vt2D; %th der. inversa
            dCcapn_n2D =  + FCapn.*(1-n2D./N2).*nF(inQW)./tauscatnT./dn2D_EFn2D./Vt2D + ...
                + (1-FCapn)./N2./tauscatnT.*nF(inQW);
        elseif(iTappo==0) % termine degli svizzeri senza il tappo
            FCapn=exp((EFn2D-EFn3D(1,:))./Vt2D);
            Ccapn = - (1-FCapn)./tauscatnT.*nF(inQW);
            dCcapn_N3D =  - (1-FCapn)./tauscatnT + ...
                - FCapn.*nF(inQW)./tauscatnT./dN3D_Efn3D./Vt2D; %th der. inversa
            dCcapn_n2D =  + FCapn.*nF(inQW)./tauscatnT./dn2D_EFn2D./Vt2D;
        elseif(iTappo==2) % termine Debernardi
            Ccapn = - ((nF(inQW)-mode.n3Di(inQW))./tauscatnT - (n2D-n2Di)./(WQW*tauescn));
            dCcapn_N3D =  - 1./tauscatnT;
            dCcapn_n2D =  + 1./(WQW*tauescn);
        elseif(iTappo==4) % termine degli svizzeri con leakage
            FCapn=exp((EFn2D-EFn3D(1,:))./Vt2D);
            Ccapn = - (1-FCapn).*(1-n2D./N2)./tauscatnT.*nF(inQW) + (n2D-n2Di)./(WQW*tauleakn);
            dCcapn_N3D =  - (1-FCapn).*(1-n2D./N2)./tauscatnT + ...
                - FCapn.*(1-n2D./N2).*nF(inQW)./tauscatnT./dN3D_Efn3D./Vt2D; %th der. inversa
            dCcapn_n2D =  + FCapn.*(1-n2D./N2).*nF(inQW)./tauscatnT./dn2D_EFn2D./Vt2D + ...
                + (1-FCapn)./N2./tauscatnT.*nF(inQW) + 1./(WQW*tauleakn);
        elseif(iTappo==3) % svizzeri senza tappo e con leakage
            FCapn=exp((EFn2D-EFn3D(1,:))./Vt2D);
            Ccapn = - (1-FCapn)./tauscatnT.*nF(inQW) + (n2D-n2Di)./(WQW*tauleakn);
            dCcapn_N3D =  - (1-FCapn)./tauscatnT + ...
                - FCapn.*nF(inQW)./tauscatnT./dN3D_Efn3D./Vt2D; %th der. inversa
            dCcapn_n2D =  + FCapn.*nF(inQW)./tauscatnT./dn2D_EFn2D./Vt2D + 1./(WQW*tauleakn);
        elseif(iTappo==6) % cattura-fuga alla Cappelluti/Gioannini maniera
            taucapn = tauscatnT;
            tauescn = taucapn.*mode.rn{indQW};
            %
            Ccapn = - (1-n2D./N2).*nF(inQW)./taucapn + n2D./(WQW*tauescn);
            dCcapn_N3D =  - (1-n2D./N2)./taucapn;
            dCcapn_n2D =  + 1./N2./taucapn.*nF(inQW) + 1./(WQW*tauescn);
        elseif(iTappo==7 | iTappo==8) % tappo esponenziale
            FCapn=exp((EFn2D-EFn3D(1,:))./Vt2D);
            eD=exp(-n2D./N2);
            Ccapn = - (1-FCapn).*eD./tauscatnT.*nF(inQW);
            dCcapn_N3D =  - (1-FCapn).*eD./tauscatnT + ...
                - FCapn.*eD.*nF(inQW)./tauscatnT./dN3D_Efn3D./Vt2D; %th der. inversa
            dCcapn_n2D =  + FCapn.*eD.*nF(inQW)./tauscatnT./dn2D_EFn2D./Vt2D + ...
                + (1-FCapn).*eD./N2./tauscatnT.*nF(inQW);
        end
        
        CcapnQW=Ccapn*mode.CarrierNorm/mode.CarrierNorm2D;
        dCcapn_N3DQW=dCcapn_N3D*mode.CarrierNorm/mode.CarrierNorm2D;
        dCcapn_n2DQW=dCcapn_n2D*mode.CarrierNorm/mode.CarrierNorm2D;
        
        CcapnM1=CcapnQW;
        CcapnM2=CcapnQW;
        
        %        Xlimi=mesh.xgrid2(3)*1.5; Original code
        Xlimi=mode.rAperture*1e-4*1.5;
        fii=find(mesh.xgrid<=Xlimi);
        fie=find(mesh.xgrid>Xlimi);
        
        CcapnM1(fie)=0;
        CcapnM2(fii)=0;
        MM1 = qel.*WQW.*[Lp1.*CcapnM1(iiQW1) Lp2.*CcapnM1(iiQW2)];
        MM2 = qel.*WQW.*[Lp1.*CcapnM2(iiQW1) Lp2.*CcapnM2(iiQW2)];
        
        %             MM = qel.*WQW.*[Lp1.*Ccapn(iiQW1) Lp2.*Ccapn(iiQW2)];
        IntCcapn=IntCcapn+[sum(MM1) sum(MM2)];  % C/s
        
        mode.IntCcapn=IntCcapn;
        
        %             xdx=[0 Lp1];
        %  figure, plot(1e4*mesh.xgrid,[JpQW_x+JnQW_x; JpQW_y+JnQW_y])
        %             Jx=JpQW_x+JnQW_x;
        %             Jy=JpQW_y+JnQW_y;
        %
        %             JnQW_x=mode.Jn_x(inQW);
        %             JnQW_y=mode.Jn_y(inQW);
        %             JpQW_x=mode.Jp_x(inQW);
        %             JpQW_y=mode.Jp_y(inQW);
        
        %
        % assembling in 2D electrons equations
        MM = qel.*WQW.*[Lp1.*CcapnQW(iiQW1) Lp2.*CcapnQW(iiQW2)];
        tvet = tvet + sparse(VV+ijrQW,1,MM,neq,1);
        MM = qel.*WQW.*[Lp1.*dCcapn_N3DQW(iiQW1) Lp2.*dCcapn_N3DQW(iiQW2)];
        Jmat0 = Jmat0 + sparse(VV+ijrQW,nn+inQW(ijrQW),MM,neq,neq);
        MM = qel.*WQW.*[Lp1.*dCcapn_n2DQW(iiQW1) Lp2.*dCcapn_n2DQW(iiQW2)];
        Jmat0 = Jmat0 + sparse(VV+ijrQW,VV+ijrQW,MM,neq,neq);
        %
        % assembling in 3D electrons equations
        Ccapn3D = zeros(1,mesh.nn);
        dCcapn3D_N3D = zeros(1,mesh.nn);
        dCcapn3D_n2D = zeros(1,mesh.nn);
        %
        Ccapn3D(IQW)=Ccapn(IIQWP(IQW));
        dCcapn3D_N3D(IQW)=dCcapn_N3D(IIQWP(IQW));
        dCcapn3D_n2D(IQW)=dCcapn_n2D(IIQWP(IQW));
        %
        MM = - qel.*[SeQW1.*Ccapn3D(in1) SeQW2.*Ccapn3D(in2) SeQW3.*Ccapn3D(in3)];
        tvet = tvet + sparse(nn+ijr,1,MM(mask_ijr),neq,1);
        MM = - qel.*[SeQW1.*dCcapn3D_N3D(in1) SeQW2.*dCcapn3D_N3D(in2) SeQW3.*dCcapn3D_N3D(in3)];
        Jmat0 = Jmat0 + sparse(nn+ijr,nn+INQWP(ijr),MM(mask_ijr),neq,neq);
        MM = - qel.*[SeQW1.*dCcapn3D_n2D(in1) SeQW2.*dCcapn3D_n2D(in2) SeQW3.*dCcapn3D_n2D(in3)];
        Jmat0 = Jmat0 + sparse(nn+ijr,VV+IIQWP(ijr),MM(mask_ijr),neq,neq);
        
        
        %             MM = [SeQW1(in1) SeQW2(in2) SeQW3(in3)];
        %             MM = [SeQW1 SeQW2 SeQW3];
        %             figure(34891),hold on
        %             plot(node(1,:),node(2,:),'b*')
        %             plot(node(1,IQW),node(2,IQW),'ro')
        RQW=zeros(1,nn);
        RQW(IQW)=R(IQW);
        %            MM = qel.*[SeQW1.*(Ccapn3D(in1)+RQW(in1)) SeQW2.*(Ccapn3D(in2)+RQW(in2)) SeQW3.*(Ccapn3D(in3)+RQW(in3))];
        %             MM = qel.*[SeQW1.*(Ccapn3D(in1)) SeQW2.*(Ccapn3D(in2)) SeQW3.*(Ccapn3D(in3))];
        %             Int_Anal=pi*(mesh.xgrid(end)-mesh.xgrid(1))^2*WQW*indQW
        %            IntRecQW3Dn=IntRecQW3Dn+sum(sparse(ijr,1,MM(mask_ijr),nn,1));
        
        
        VaIn=RQW;
        MM = qel.*[SeQW1.*VaIn(in1) SeQW2.*VaIn(in2) SeQW3.*VaIn(in3)];
        IntRecQW3Dn=IntRecQW3Dn+sum(sparse(ijr,1,MM(mask_ijr),nn,1));
        
        %             ErrRelInt=abs(Int_Anal-IntegralQW)./abs(Int_Anal)
        
        mode.IntRecN=IntRecQW3Dn;
        
        dP3D_Efp3D =  - dpF(inQW)./Vt2D;
        
        if(iTappo==1) % esattamente il termine degli svizzeri
            FCapp=exp(-(EFp2D-EFp3D(1,:))./Vt2D);
            Ccapp = - (1-FCapp).*(1-p2D./P2)./tauscatpT.*pF(inQW);
            dCcapp_P3D =  - (1-FCapp).*(1-p2D./P2)./tauscatpT + ...
                + FCapp.*(1-p2D./P2).*pF(inQW)./tauscatpT./dP3D_Efp3D./Vt2D; %th der. inversa
            dCcapp_p2D =  - FCapp.*(1-p2D./P2).*pF(inQW)./tauscatpT./dp2D_EFp2D./Vt2D + ...
                + (1-FCapp)./P2./tauscatpT.*pF(inQW);
        elseif(iTappo==0) % termine degli svizzeri senza il tappo
            FCapp=exp(-(EFp2D-EFp3D(1,:))./Vt2D);
            Ccapp = - (1-FCapp)./tauscatpT.*pF(inQW);
            dCcapp_P3D =  - (1-FCapp)./tauscatpT + ...
                + FCapp.*pF(inQW)./tauscatpT./dP3D_Efp3D./Vt2D; %th der. inversa
            dCcapp_p2D =  - FCapp.*pF(inQW)./tauscatpT./dp2D_EFp2D./Vt2D;
        elseif(iTappo==2) % termine Debernardi
            Ccapp = - ((pF(inQW)-mode.p3Di(inQW))./tauscatpT - (p2D-p2Di)./(WQW*tauescp));
            dCcapp_P3D =  - 1./tauscatpT;
            dCcapp_p2D =  + 1./(WQW*tauescp);
        elseif(iTappo==4) % termine degli svizzeri con leakage
            FCapp=exp(-(EFp2D-EFp3D(1,:))./Vt2D);
            Ccapp = - (1-FCapp).*(1-p2D./P2)./tauscatpT.*pF(inQW) + (p2D-p2Di)./(WQW*tauleakp);
            dCcapp_P3D =  - (1-FCapp).*(1-p2D./P2)./tauscatpT + ...
                + FCapp.*(1-p2D./P2).*pF(inQW)./tauscatpT./dP3D_Efp3D./Vt2D; %th der. inversa
            dCcapp_p2D =  - FCapp.*(1-p2D./P2).*pF(inQW)./tauscatpT./dp2D_EFp2D./Vt2D + ...
                + (1-FCapp)./P2./tauscatpT.*pF(inQW) + 1./(WQW*tauleakp);
        elseif(iTappo==3) % svizzeri senza tappo e con leakage
            FCapp=exp(-(EFp2D-EFp3D(1,:))./Vt2D);
            Ccapp = - (1-FCapp)./tauscatpT.*pF(inQW) + (p2D-p2Di)./(WQW*tauleakp);
            dCcapp_P3D =  - (1-FCapp)./tauscatpT + ...
                + FCapp.*pF(inQW)./tauscatpT./dP3D_Efp3D./Vt2D; %th der. inversa
            dCcapp_p2D =  - FCapp.*pF(inQW)./tauscatpT./dp2D_EFp2D./Vt2D + 1./(WQW*tauleakp);
        elseif(iTappo==6) % cattura-fuga alla Cappelluti/Gioannini maniera
            taucapp = tauscatpT;
            tauescp = taucapp.*mode.rp{indQW};
            %
            Ccapp = - (1-p2D./P2).*pF(inQW)./taucapp + p2D./(WQW.*tauescp);
            dCcapp_P3D =  - (1-p2D./P2)./taucapp;
            dCcapp_p2D =  + 1./P2./taucapp.*pF(inQW) + 1./(WQW*tauescp);
        elseif(iTappo==7 | iTappo==8) % mio
            pD=exp(-p2D./P2);
            FCapp=exp(-(EFp2D-EFp3D(1,:))./Vt2D);
            Ccapp = - (1-FCapp).*pD./tauscatpT.*pF(inQW);
            dCcapp_P3D =  - (1-FCapp).*pD./tauscatpT + ...
                + FCapp.*pD.*pF(inQW)./tauscatpT./dP3D_Efp3D./Vt2D; %th der. inversa
            dCcapp_p2D =  - FCapp.*pD.*pF(inQW)./tauscatpT./dp2D_EFp2D./Vt2D + ...
                + (1-FCapp).*pD./P2./tauscatpT.*pF(inQW);
            
        end
        
        CcappQW=Ccapp*mode.CarrierNorm/mode.CarrierNorm2D;
        dCcapp_P3DQW=dCcapp_P3D*mode.CarrierNorm/mode.CarrierNorm2D;
        dCcapp_p2DQW=dCcapp_p2D*mode.CarrierNorm/mode.CarrierNorm2D;
        
        CcappM1=CcappQW;
        CcappM2=CcappQW;
        CcappM1(fie)=0;
        CcappM2(fii)=0;
        MM1 = qel.*WQW.*[Lp1.*CcappM1(iiQW1) Lp2.*CcappM1(iiQW2)];
        MM2 = qel.*WQW.*[Lp1.*CcappM2(iiQW1) Lp2.*CcappM2(iiQW2)];
        
        IntCcapp=IntCcapp+[sum(MM1) sum(MM2)];  % C/s
        mode.IntCcapp=IntCcapp;
        %
        % assembling in 2D holes equations
        MM = qel.*WQW.*[Lp1.*CcappQW(iiQW1) Lp2.*CcappQW(iiQW2)];
        tvet = tvet + sparse(WW+ijrQW,1,MM,neq,1);
        MM = qel.*WQW.*[Lp1.*dCcapp_P3DQW(iiQW1) Lp2.*dCcapp_P3DQW(iiQW2)];
        Jmat0 = Jmat0 + sparse(WW+ijrQW,pp+inQW(ijrQW),MM,neq,neq);
        MM = qel.*WQW.*[Lp1.*dCcapp_p2DQW(iiQW1) Lp2.*dCcapp_p2DQW(iiQW2)];
        Jmat0 = Jmat0 + sparse(WW+ijrQW,WW+ijrQW,MM,neq,neq);
        %
        % assembling in 3D holes equations
        Ccapp3D = zeros(1,mesh.nn);
        dCcapp3D_N3D = zeros(1,mesh.nn);
        dCcapp3D_n2D = zeros(1,mesh.nn);
        %
        Ccapp3D(IQW)=Ccapp(IIQWP(IQW));
        dCcapp3D_N3D(IQW)=dCcapp_P3D(IIQWP(IQW));
        dCcapp3D_n2D(IQW)=dCcapp_p2D(IIQWP(IQW));
        %
        MM = - qel.*[SeQW1.*Ccapp3D(in1) SeQW2.*Ccapp3D(in2) SeQW3.*Ccapp3D(in3)];
        tvet = tvet + sparse(pp+ijr,1,MM(mask_ijr),neq,1);
        MM = - qel.*[SeQW1.*dCcapp3D_N3D(in1) SeQW2.*dCcapp3D_N3D(in2) SeQW3.*dCcapp3D_N3D(in3)];
        Jmat0 = Jmat0 + sparse(pp+ijr,pp+INQWP(ijr),MM(mask_ijr),neq,neq);
        MM = - qel.*[SeQW1.*dCcapp3D_n2D(in1) SeQW2.*dCcapp3D_n2D(in2) SeQW3.*dCcapp3D_n2D(in3)];
        Jmat0 = Jmat0 + sparse(pp+ijr,WW+IIQWP(ijr),MM(mask_ijr),neq,neq);
        
        %            VaIn=Ccapp3D+RQW;
        %            MM = qel.*[SeQW1.*(Ccapp3D(in1)+RQW(in1)) SeQW2.*(Ccapp3D(in2)+RQW(in2)) SeQW3.*(Ccapp3D(in3)+RQW(in3))];
        %            MM = qel.*[SeQW1.*VaIn(in1) SeQW2.*VaIn(in2) SeQW3.*VaIn(in3)];
        
        VaIn=RQW;
        MM = qel.*[SeQW1.*VaIn(in1) SeQW2.*VaIn(in2) SeQW3.*VaIn(in3)];
        
        IntRecQW3Dp=IntRecQW3Dp+sum(sparse(ijr,1,MM(mask_ijr),nn,1));
        
        mode.IntRecP=IntRecQW3Dp;
        
        %
        % Recombination of bound carriers
        R2D = 0; dR2D_n2D=0; dR2D_p2D=0;
        np2D=n2D.*p2D-n2Di.*p2Di;
        %
        % Shockley-Read-Hall recombination term
        taunSRH=mesh.taunQW(inQW);
        taupSRH=mesh.taupQW(inQW);
        denSRH2D = taupSRH.*(n2D+n2Di)+taunSRH.*(p2D+p2Di);
        RSRH2D = np2D./denSRH2D;
        dRSRH_n2D = (p2D.*denSRH2D - taupSRH.*np2D)./(denSRH2D.^2);
        dRSRH_p2D = (n2D.*denSRH2D - taunSRH.*np2D)./(denSRH2D.^2);
        R2D = R2D + RSRH2D;
        dR2D_n2D = dR2D_n2D + dRSRH_n2D;
        dR2D_p2D = dR2D_p2D + dRSRH_p2D;
        
        VI=RSRH2D;
        MM = qel.*[Lp1.*VI(iiQW1) Lp2.*VI(iiQW2)];
        IntSRH=IntSRH+sum(MM);
        
        %
        % Radiative recombination term
        % Semiconductor Bloch Equations (SBE) model
        %            nQW=n2D./WQW;
        %            pQW=p2D./WQW;
        %            np=(nQW+pQW)/2-(n2Di+p2Di)/2/WQW;
        %            IntRecQW3Dn
        
        np_E=n2D-n2Di;
        np_H=p2D-p2Di;
        
        np_E=n2D;
        np_H=p2D;
        
        [Rsp,dRspE,dRspH]=f_InterpRsp_lin(np_E,np_H,indQW);
        
        Rsp=Rsp*mode.fat_gain;
        dRspE=dRspE*mode.fat_gain;
        dRspH=dRspH*mode.fat_gain;
        %             if mode.CarrierNorm>1
        %                 Rrad2D_SBE = Rsp; % from 1/(s*cm^3) to 1/(s*cm^2)
        %                 Rrad2D = Rsp; % from 1/(s*cm^3) to 1/(s*cm^2)
        %                 dRrad2D_n2D = dRspE; % 1/cm and cm simplify
        %                 dRrad2D_p2D = dRspH; % 1/cm and cm simplify
        %             else
        Rrad2D_SBE = Rsp*WQW; % from 1/(s*cm^3) to 1/(s*cm^2)
        Rrad2D = Rsp*WQW; % from 1/(s*cm^3) to 1/(s*cm^2)
        dRrad2D_n2D = dRspE*WQW; % 1/cm and cm simplify
        dRrad2D_p2D = dRspH*WQW; % 1/cm and cm simplify
        %             end
        
        VI=Rrad2D;
        MM = qel.*[Lp1.*VI(iiQW1) Lp2.*VI(iiQW2)];
        IntRad=IntRad+sum(MM);
        mode.IRsp(mode.indv,indQW) = sum (MM);
        
        MM = [Lp1.*Rrad2D(iiQW1) Lp2.*Rrad2D(iiQW2)];
        frsp=mode.frsp;
        Psp=Psp+frsp*1000*sum(sparse(MM))*h*(Clight*1e-2/mean(1e-9*mode.vlambda)); % milliwatt
        
        % "Brad" model
        %         Brad2D = mesh.brad(inQW)./WQW; % Brad in cm^2/s (/WQW!)
        %         Rrad2D = Brad2D.*np2D; % cm^(-2)/s
        %         dRrad2D_n2D = Brad2D.*p2D;
        %         dRrad2D_p2D = Brad2D.*n2D;
        %
        R2D = R2D + Rrad2D;
        dR2D_n2D = dR2D_n2D + dRrad2D_n2D;
        dR2D_p2D = dR2D_p2D + dRrad2D_p2D;
        %
        % Auger recombination model
        %             if mode.CarrierNorm>1
        %                 Cnnp2D=mode.FatAuger23D*mesh.Cnnp(inQW);
        %                 Cppn2D=mode.FatAuger23D*mesh.Cppn(inQW);
        %             else
        Cnnp2D=mode.FatAuger23D*mesh.Cnnp(inQW)./(WQW.^2);
        Cppn2D=mode.FatAuger23D*mesh.Cppn(inQW)./(WQW.^2);
        %             end
        RAuger_n2D = Cnnp2D.*n2D.*np2D;
        RAuger_p2D = Cppn2D.*p2D.*np2D;
        dRAuger_n2D_n2D = 2*Cnnp2D.*np2D; % derivative of RAuger_n2D w.r.t. n2D
        dRAuger_n2D_p2D = Cnnp2D.*n2D.*n2D;  % derivative of RAuger_n2D w.r.t. p2D
        dRAuger_p2D_n2D = Cppn2D.*p2D.*p2D; % derivative of RAuger_p2D w.r.t. n2D
        dRAuger_p2D_p2D = 2*Cppn2D.*np2D;  % derivative of RAuger_p2D w.r.t. p2D
        RAuger2D = RAuger_n2D + RAuger_p2D;
        dRAuger_n2D = dRAuger_n2D_n2D + dRAuger_p2D_n2D;
        dRAuger_p2D = dRAuger_n2D_p2D + dRAuger_p2D_p2D;
        
        %             dRAuger_n2D = 2*Cnnp2D.*np2D + Cppn2D.*p2D.*p2D;
        %             dRAuger_p2D = 2*Cppn2D.*np2D + Cnnp2D.*n2D.*n2D;
        
        VI=RAuger2D;
        MM = qel.*[Lp1.*VI(iiQW1) Lp2.*VI(iiQW2)];
        IntAug=IntAug+sum(MM);
        
        R2D = R2D + RAuger2D;
        dR2D_n2D = dR2D_n2D + dRAuger_n2D;
        dR2D_p2D = dR2D_p2D + dRAuger_p2D;
        
        % Auger expulsion model, from 2013Deppner_SPIE
        if(mode.AugerExpulsion==1)
            % Assembling n2D expulsion term (n2D equation only)
            MM = qel.*[Lp1.*RAuger_n2D(iiQW1) Lp2.*RAuger_n2D(iiQW2)];
            tvet = tvet + sparse(VV+ijrQW,1,MM,neq,1);
            % Assembling n2D expulsion jacobian w.r.t. n2D
            MM = qel.*[Lp1.*dRAuger_n2D_n2D(iiQW1) Lp2.*dRAuger_n2D_n2D(iiQW2)];
            Jmat0 = Jmat0 + sparse(VV+ijrQW,VV+ijrQW,MM,neq,neq);
            % Assembling n2D expulsion jacobian w.r.t. p2D
            MM = qel.*[Lp1.*dRAuger_n2D_p2D(iiQW1) Lp2.*dRAuger_n2D_p2D(iiQW2)];
            Jmat0 = Jmat0 + sparse(VV+ijrQW,WW+ijrQW,MM,neq,neq);
            %
            % Assembling p2D expulsion term (p2D equation only)
            MM = qel.*[Lp1.*RAuger_p2D(iiQW1) Lp2.*RAuger_p2D(iiQW2)];
            tvet = tvet + sparse(WW+ijrQW,1,MM,neq,1);
            % Assembling p2D expulsion jacobian w.r.t. p2D
            MM = qel.*[Lp1.*dRAuger_p2D_p2D(iiQW1) Lp2.*dRAuger_p2D_p2D(iiQW2)];
            Jmat0 = Jmat0 + sparse(WW+ijrQW,WW+ijrQW,MM,neq,neq);
            % Assembling p2D expulsion jacobian w.r.t. n2D
            MM = qel.*[Lp1.*dRAuger_p2D_n2D(iiQW1) Lp2.*dRAuger_p2D_n2D(iiQW2)];
            Jmat0 = Jmat0 + sparse(WW+ijrQW,VV+ijrQW,MM,neq,neq);
            %
            % Auger regeneration: Gaussian profile
            Deltaz_E=190e-7;
            Deltaz_H=170e-7;
            
            B=mode.Auger_broad;
            f_reg=mode.Fat_regeneration;
            
            GaussEnvelopE=exp(-2*((mesh.ygrid-(yQW+Deltaz_E))./(B*1e-7)).^2);
            GaussEnvelopH=exp(-2*((mesh.ygrid-(yQW-Deltaz_H))./(B*1e-7)).^2);
            GaussEnvelop=exp(-2*((mesh.ygrid-yQW)./(B*1e-7)).^2);
            GaSum=GaSum+GaussEnvelop;
            %                pausak
            % Norm2=sqrt(pi/2)*B*1e-7;
            Norm2=trapz(mesh.ygrid,GaussEnvelop);
            Norm2E=trapz(mesh.ygrid,GaussEnvelopE);
            Norm2H=trapz(mesh.ygrid,GaussEnvelopH);
            GaussEnvelop=GaussEnvelop./Norm2;
            GaussEnvelopE=GaussEnvelopE/Norm2E;
            GaussEnvelopH=GaussEnvelopH/Norm2H;
            %
            % Assembling n2D regeneration term (n3D equation only)
            GAuger_n2D      = - f_reg*RAuger_n2D;
            GAuger_n3D      = reshape(GaussEnvelopE.'*GAuger_n2D,1,nn);
            MM              = qel.*[Se1.*GAuger_n3D(in1) Se2.*GAuger_n3D(in2) Se3.*GAuger_n3D(in3)];
            tvet            = tvet + sparse(nn+ijs,1,MM(mask_ijs),neq,1);
            % Assembling n3D regeneration jacobian w.r.t. n2D
            dGAuger_n2D_n2D = - f_reg*dRAuger_n2D_n2D;
            dGAuger_n3D_n2D = reshape(GaussEnvelopE.'*dGAuger_n2D_n2D,1,nn);
            MM              = qel.*[Se1.*dGAuger_n3D_n2D(in1) Se2.*dGAuger_n3D_n2D(in2) Se3.*dGAuger_n3D_n2D(in3)];
            Jmat0           = Jmat0 + sparse(nn+ijs,VV+IIQWP(ijs),MM(mask_ijs),neq,neq);
            % Assembling n3D regeneration jacobian w.r.t. p2D
            dGAuger_n2D_p2D = - f_reg*dRAuger_n2D_p2D;
            dGAuger_n3D_p2D = reshape(GaussEnvelopE.'*dGAuger_n2D_p2D,1,nn);
            MM              = qel.*[Se1.*dGAuger_n3D_p2D(in1) Se2.*dGAuger_n3D_p2D(in2) Se3.*dGAuger_n3D_p2D(in3)];
            Jmat0           = Jmat0 + sparse(nn+ijs,WW+IIQWP(ijs),MM(mask_ijs),neq,neq);
            %
            % Assembling p2D regeneration term (p3D equation only)
            GAuger_p2D      = - f_reg*RAuger_p2D;
            GAuger_p3D      = reshape(GaussEnvelopH.'*GAuger_p2D,1,nn);
            MM              = qel.*[Se1.*GAuger_p3D(in1) Se2.*GAuger_p3D(in2) Se3.*GAuger_p3D(in3)];
            tvet            = tvet + sparse(pp+ijs,1,MM(mask_ijs),neq,1);
            % Assembling p3D regeneration jacobian w.r.t. p2D
            dGAuger_p2D_p2D = - f_reg*dRAuger_p2D_p2D;
            dGAuger_p3D_p2D = reshape(GaussEnvelopH.'*dGAuger_p2D_p2D,1,nn);
            MM              = qel.*[Se1.*dGAuger_p3D_p2D(in1) Se2.*dGAuger_p3D_p2D(in2) Se3.*dGAuger_p3D_p2D(in3)];
            Jmat0           = Jmat0 + sparse(pp+ijs,WW+IIQWP(ijs),MM(mask_ijs),neq,neq);
            % Assembling p3D regeneration jacobian w.r.t. n2D
            dGAuger_p2D_n2D = - f_reg*dRAuger_p2D_n2D;
            dGAuger_p3D_n2D = reshape(GaussEnvelopH.'*dGAuger_p2D_n2D,1,nn);
            MM              = qel.*[Se1.*dGAuger_p3D_n2D(in1) Se2.*dGAuger_p3D_n2D(in2) Se3.*dGAuger_p3D_n2D(in3)];
            Jmat0           = Jmat0 + sparse(pp+ijs,VV+IIQWP(ijs),MM(mask_ijs),neq,neq);
            
        end
        
        %
        % Additional SRH-like leakage terms
        if(mode.leakage>0)
            % tauleakn
            % tauleakp
            %                tauleakn=1e-8;
            %                tauleakp=1e-8;
            tauleakn=Fat_capleak;
            tauleakp=Fat_capleak;
            %                 RLeakage2D=(n2D-n2Di)./tauleakn + (p2D-p2Di)./tauleakp;
            %                 dRLeakage_n2D = 1./tauleakn;
            %                 dRLeakage_p2D = 1./tauleakp;
            %            np2D=n2D.*p2D-n2Di.*p2Di;
            denLeakage2D = tauleakp.*(n2D+n2Di)+tauleakn.*(p2D+p2Di);
            RLeakage2D = np2D./denLeakage2D;
            VI=RLeakage2D;
            MM = qel.*[Lp1.*VI(iiQW1) Lp2.*VI(iiQW2)];
            IntLea=IntLea+sum(MM);
            
            dRLeakage_n2D = (p2D.*denLeakage2D - tauleakp.*np2D)./(denLeakage2D.^2);
            dRLeakage_p2D = (n2D.*denLeakage2D - tauleakn.*np2D)./(denLeakage2D.^2);
            
            R2D = R2D + RLeakage2D;
            dR2D_n2D = dR2D_n2D + dRLeakage_n2D;
            dR2D_p2D = dR2D_p2D + dRLeakage_p2D;
        else
            RLeakage2D=zeros(1,mesh.nnx);
        end
        %
        % 2D electron equation recombination term ---------------------------------
        MM = qel.*[Lp1.*R2D(iiQW1) Lp2.*R2D(iiQW2)];
        tvet = tvet + sparse(VV+ijrQW,1,MM,neq,1);
        mode.IR2D(mode.indv,indQW)=sum(MM) ;
        MM = qel.*[Lp1.*dR2D_p2D(iiQW1) Lp2.*dR2D_p2D(iiQW2)];
        Jmat0 = Jmat0 + sparse(VV+ijrQW,WW+ijrQW,MM,neq,neq);
        MM = qel.*[Lp1.*dR2D_n2D(iiQW1) Lp2.*dR2D_n2D(iiQW2)];
        Jmat0 = Jmat0 + sparse(VV+ijrQW,VV+ijrQW,MM,neq,neq);
        %
        % 2D hole equation recombination term -------------------------------------
        MM = qel.*[Lp1.*R2D(iiQW1) Lp2.*R2D(iiQW2)];
        tvet = tvet + sparse(WW+ijrQW,1,MM,neq,1);
        MM = qel.*[Lp1.*dR2D_p2D(iiQW1) Lp2.*dR2D_p2D(iiQW2)];
        Jmat0 = Jmat0 + sparse(WW+ijrQW,WW+ijrQW,MM,neq,neq);
        MM = qel.*[Lp1.*dR2D_n2D(iiQW1) Lp2.*dR2D_n2D(iiQW2)];
        Jmat0 = Jmat0 + sparse(WW+ijrQW,VV+ijrQW,MM,neq,neq);
        %
        mode.RSRHQW(IQW)=mode.RSRHQW(IQW)+RSRH2D(IIQWP(IQW))./WQW;
        mode.RradQW(IQW)=mode.RradQW(IQW)+Rrad2D(IIQWP(IQW))./WQW;
        mode.RradQW_SBE(IQW)=mode.RradQW_SBE(IQW)+Rrad2D_SBE(IIQWP(IQW))./WQW;
        mode.RAugerQW(IQW)=mode.RAugerQW(IQW)+RAuger2D(IIQWP(IQW))./WQW;
        mode.RLeakageQW(IQW)=mode.RLeakageQW(IQW)+RLeakage2D(IIQWP(IQW))./WQW;
        %
        %x=mesh.xgrid;
        %		Cn=sum(reshape(mode.Ccapn3D,mesh.nny,mesh.nnx));
        %		Cp=sum(reshape(mode.Ccapp3D,mesh.nny,mesh.nnx));
        % Icn=qel*2*pi*trapz(x,x.*Cn)*1000*WQW;
        % Icp=qel*2*pi*trapz(x,x.*Cp)*1000*WQW;
        %
        mode.Em{indQW}=Em; mode.n2D{indQW}=n2D; mode.EFn2D{indQW,mode.ind_v0}=EFn2D;
        mode.Ed{indQW}=Ed; mode.p2D{indQW}=p2D; mode.EFp2D{indQW,mode.ind_v0}=EFp2D;
        mode.N2{indQW}=N2; mode.P2{indQW}=P2;
        mode.Ccapn{indQW,mode.ind_v0}=Ccapn; mode.Ccapp{indQW,mode.ind_v0}=Ccapp;
        %
        mode.Ccapn3D=mode.Ccapn3D-Ccapn3D; mode.Ccapp3D=mode.Ccapp3D-Ccapp3D;
        %
        % =========================================================================
        % Assembling photon rate equations: stimulated emission
        % =========================================================================
        vph=Clight./mode.nindexQW; % phase velocity in GaAs quantum well
        TQW=mesh.T(inQW);
        diffMQW=find(abs(diff(cell2mat(mesh.MQWcell)))>2);
        %
        for indMode=1:nmodes
            Gamma_z=mode.Gamma_z(indMode,indQW);
            if isempty(diffMQW)
                % 1 MQW region
                Campo_attivo=mode.Gamma_z(indMode,:)./mode.Gamma_z(indMode,ceil(end/2));
            else
                % more than one MQW region
                Campo_attivo=mode.Gamma_z(indMode,:)./mode.Gamma_z(indMode,ceil(diffMQW(1)/2));
            end
            Fattore_attivo=Campo_attivo(indQW);
            
            [g,dgE,dgH,rsp,drspE,drspH] = f_InterpGain_lin(n2D,p2D,indQW,indMode);
            
            %
            if isfield(mode,'fat_gainG')
                Fat_soloGain=mode.fat_gainG;
            else
                Fat_soloGain=1;
            end
            
            if isfield(mode,'nlG')
                nlG=mode.nlG;
            else
                nlG=1;
            end
            
            g=vph*mode.fat_gain*Fat_soloGain*g.*nlG;
            
            mode.g{indQW,indMode}=g; % saving just for the last mode..
            
            dgE=vph*mode.fat_gain*Fat_soloGain*dgE.*nlG;
            dgH=vph*mode.fat_gain*Fat_soloGain*dgH.*nlG;
            
            rsp=rsp*mode.fat_gain;
            drspE=drspE*mode.fat_gain;
            drspH=drspH*mode.fat_gain;
            %
            %-- Stimulated emission: rate equation
            mode.matgain=mode.matgain+g/vph; % saving material gain for VELM
            %            mode.DeltaN=mode.DeltaN+DeltaN;
            %
            E2 = mode.E2(indMode,1:nnQW); % electric field intensity (normalized)
            %
            FLos=mode.FLos;
            Lm = FLos*mode.Lm(indMode); % losses, 1/s
            %
            gE = g.*E2; % gain-field product
            MM = [Lp1.*gE(iiQW1) Lp2.*gE(iiQW2)];
            
            mode.gE(mode.indv,indQW,indMode,:)=gE;
            
            gm = sum(diag(sparse(ijrQW,ijrQW,MM,nnQW,nnQW)));
            dgm_nE = dgE.*E2;
            dgm_nH = dgH.*E2;
            %
            rspE = rsp.*E2;
            MM = [Lp1.*rspE(iiQW1) Lp2.*rspE(iiQW2)];
            rspMod = sum(diag(sparse(ijrQW,ijrQW,MM,nnQW,nnQW)));
            dRsp_nE = drspE.*E2;
            dRsp_nH = drspH.*E2;
            %
            % Assembling spontaneous emission power
            Sq=Gamma_z.*mode.fPES(indMode).*rspMod;
            mode.Sq(mode.ind_v0)=Sq;
            tvet(ss+indMode)=tvet(ss+indMode)+Sq;
            %
            % Assembling derivatives of photon equation w.r.t. carriers
            PJacobE = Gamma_z.*dgm_nE*Pst(indMode) + dRsp_nE.*Gamma_z.*mode.fPES(indMode);
            MM = [Lp1.*PJacobE(iiQW1) Lp2.*PJacobE(iiQW2)];
            Jmat0 = Jmat0+sparse(ss+indMode,VV+ijrQW,MM,neq,neq);
            PJacobH = Gamma_z.*dgm_nH*Pst(indMode) + dRsp_nH.*Gamma_z.*mode.fPES(indMode);
            MM = [Lp1.*PJacobH(iiQW1) Lp2.*PJacobH(iiQW2)];
            Jmat0 = Jmat0+sparse(ss+indMode,WW+ijrQW,MM,neq,neq);
            %
            % Assembling derivative w.r.t. Pst of photon rate equation (linear)
            MM = (Gamma_z.*gm-Lm/NQW);    % la somma su tutti i QW deve dare Lm
            MMgain=MM+MMgain;
            Kmat0=Kmat0+sparse(ss+indMode,ss+indMode,MM,neq,neq);
            %
            % for display only
            mode.gmod(indMode)=mode.gmod(indMode)+Gamma_z.*gm./vph;
            mode.lmod(indMode)=Lm./vph;
            %
            % Equazione Portatori 2D
            %-- Stimulated emission: recombination term
            Rst = gE.*Pst(indMode).*Fattore_attivo*mode.fPdif(indMode)/mode.CarrierNorm2D*WQW;
            % 2D electron equation recombination term ---------------------------------
            MM = qel.*[Lp1.*Rst(iiQW1) Lp2.*Rst(iiQW2)];
            mode.Irst(mode.indv,indQW)=sum(MM) ;
            IntStim=IntStim+sum(MM);
            
            tvet = tvet + sparse(VV+ijrQW,1,MM,neq,1);
            tvet = tvet + sparse(WW+ijrQW,1,MM,neq,1);
            %
            % Derivative of electrons equation Rst with respect to electrons
            dRst_n2D = Fattore_attivo*dgm_nE.*Pst(indMode).*mode.fPdif(indMode)/mode.CarrierNorm2D*WQW;
            MM = qel.*[Lp1.*dRst_n2D(iiQW1) Lp2.*dRst_n2D(iiQW2)];
            Jmat0 = Jmat0+sparse(VV+ijrQW,VV+ijrQW,MM,neq,neq);
            %
            % Derivative of holes equation Rst with respect to holes
            dRst_p2D = Fattore_attivo*dgm_nH.*Pst(indMode).*mode.fPdif(indMode)/mode.CarrierNorm2D*WQW;
            MM = qel.*[Lp1.*dRst_p2D(iiQW1) Lp2.*dRst_p2D(iiQW2)];
            Jmat0 = Jmat0+sparse(WW+ijrQW,WW+ijrQW,MM,neq,neq);
            %
            % Derivative of electrons equation Rst with respect to holes
            MM = qel.*[Lp1.*dRst_p2D(iiQW1) Lp2.*dRst_p2D(iiQW2)];
            Jmat0 = Jmat0+sparse(VV+ijrQW,WW+ijrQW,MM,neq,neq);
            %
            % Derivative of holes equation with Rst respect to electrons
            MM = qel.*[Lp1.*dRst_n2D(iiQW1) Lp2.*dRst_n2D(iiQW2)];
            Jmat0 = Jmat0+sparse(WW+ijrQW,VV+ijrQW,MM,neq,neq);
            %
            % Derivatives of electrons and holes Rst equation w.r.t. Pst
            %Rst = gE.*Pst(indMode).*mode.fPdif(indMode)*WQW;
            dRst_Pi = Fattore_attivo*gE.*mode.fPdif(indMode)*WQW/mode.CarrierNorm2D;
            MM = qel.*[Lp1.*dRst_Pi(iiQW1) Lp2.*dRst_Pi(iiQW2)];
            Jmat0 = Jmat0 + sparse(VV+ijrQW,ss+indMode,MM,neq,neq);
            Jmat0 = Jmat0 + sparse(WW+ijrQW,ss+indMode,MM,neq,neq);
            
        end   %modes
        
        tvetQW(VV+[1:mesh.nnxQW{indQW}]) = tvet(VV+[1:mesh.nnxQW{indQW}]);
        tvetQW(WW+[1:mesh.nnxQW{indQW}]) = tvet(WW+[1:mesh.nnxQW{indQW}]);
        
    end  % MQW
    
    mode.MMgain(mode.ind_v0)=MMgain;
    
    curQW = ImatQW*uvet'+tvetQW;
    
    iiQW=1:nnQW;
    for indQW=1:NQW
        WQW = mesh.vWMQW{indQW}; %-- cm
        VV = vv+(indQW-1)*length(mesh.xgrid);
        WW = ww+(indQW-1)*length(mesh.xgrid);
        mode.JnQW{indQW} = curQW(VV+iiQW)/WQW;
        mode.JpQW{indQW} = curQW(WW+iiQW)/WQW;
        
    end
    %         figure,plot(mesh.xgrid,mode.JnQW{1},mesh.xgrid,mode.Jn_x(mesh.inMQW{indQW}),'--'),title('Electron currents'),legend('2D','3D')
    %         figure,plot(mesh.xgrid,mode.JpQW{1},mesh.xgrid,mode.Jp_x(mesh.inMQW{indQW}),'--'),title('Hole currents'),legend('2D','3D')
    Scheck=mode.gmod'./mode.lmod; % stimulated emission check
    %
    % Debug prints
    for indMode=1:nmodes
        fprintf('Pst = %.5e\n',Pst(indMode).')
        fprintf('Gamma_z*Gm/Lm = %.8e\n\n',full(Scheck(indMode)))
    end
    
    mode.N2D=N2Dtot; mode.P2D=P2Dtot;
    %
    % saving variables
    mode.Pst=Pst; mode.Scheck=Scheck;
    %mode.nQW{mode.ind_v0}=nQW; mode.pQW{mode.ind_v0}=pQW;
    mode.matgain=mode.matgain/(NQW*nmodes);
    mode.Psp=Psp;
    
    % assembling additional equations for passivation
    if(mesh.nnx-nnQW>0)
        iiP = nnQW+1:mesh.nnx;
        for indQW=1:NQW
            VV = vv+(indQW-1)*length(mesh.xgrid);
            WW = ww+(indQW-1)*length(mesh.xgrid);
            Kmat0 = Kmat0 + sparse(VV+iiP,VV+iiP,1,neq,neq);
            Kmat0 = Kmat0 + sparse(WW+iiP,WW+iiP,1,neq,neq);
        end
    end
    
end
%mode.DeltaN=mode.DeltaN/(NQW*nmodes);
%
% =============================================================================================100
% Apply boundary conditions, Poisson equation: phi-V-phi_r=0
% =============================================================================================100
ii=find(not(maskr));
Kmat0=Kmat0+sparse(ii,ii,1,neq,neq); % phi
%
for ic=1:nm
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
        otherwise, error('contact type unknown!')
    end
end

tvet3=sparse(neq,1);
Jmat3=sparse(neq,neq);

if mode.Shunt==1
    
    Ncol=mesh.nnx ; % 2
    mask_shunt_r=zeros(1,mesh.nny);
    Ncol=min(Ncol,10) ;
    mask_shunt_r(2:mesh.nny-1)=1;
    mask_shunt=zeros(1,nn);
    nnn=mesh.nny;
    
    indAperture=find(mesh.xgrid*1e4<=mode.rAperture*1.01);
    indAperture=mesh.nnx;
    
    for indCol=1:length(indAperture)
        mask_shunt(1+nnn*(indCol-1):nnn*indCol)=mask_shunt_r;
    end
    
    %
    Rshunt=zeros(1,nn);
    dRshunt=zeros(1,nn);
    mask_shunt=logical(mask_shunt);
    Rshunt(mask_shunt(iq))=mode.K.*(phi(mask_shunt(iq))-mode.phi_0(mask_shunt(iq)))/mode.CarrierNorm;
    dRshunt(mask_shunt(iq))=mode.K/mode.CarrierNorm;
    
    MM=qel.*[Rshunt(in1).*Se1 Rshunt(in2).*Se2 Rshunt(in3).*Se3];
    tvet3=tvet3+sparse(nn+ijs,1,MM(mask_ijs),neq,1);
    
    MM=qel.*[dRshunt(in1()).*Se1 dRshunt(in2).*Se2 dRshunt(in3).*Se3];
    Jmat3=Jmat3+sparse(nn+ijs,ijs,MM(mask_ijs),neq,neq);
    
end


% =============================================================================================100
% Apply boundary conditions, continuity equations
% =============================================================================================100
ii=find(not(masks));
if(mode.nflg), Jmat0=Jmat0+sparse(ii+nn,ii+nn,1,neq,neq); end
if(mode.pflg), Jmat0=Jmat0+sparse(ii+pp,ii+pp,1,neq,neq); end
%
% =============================================================================================100
% Compute residual
% =============================================================================================100

rvet = (Kmat0 + Jmat2) * uvet.' + tvet + tvet3;
Jmat0 = Kmat0 + Jmat0 + Jmat2 + Jmat3;


rvetTot=rvet+tvet3;
JmatTot=Jmat0+Jmat3;

if mode.oflg
    Intot=sum([IntAug IntSRH IntRad IntStim IntLea]);
    IntCa=sum(IntCcapp);%/WQW
    
    fprintf('           GR processes integrals (QW)\n')
    fprintf('    QW recombinations     |             Capture             \n')
    fprintf(['Int. Tot. = ',num2str(Intot), '    |   Int. Cap. = ',num2str(IntCa),'\n\n'])
    
    mode.IntRec=[IntSRH IntRad IntAug IntStim IntLea];
end
if(mode.AugerExpulsion==1)
    mode.GaSum=GaSum;
    %                  figure(489330), plot(mesh.ygrid*1e4,GaSum),
    %                  xlim([354 356])
end
