
function mesh=rectmesh(geom)

%
% History: created Francesco Bertazzi, still in writing!
% Description: generates a rectangular mesh as the tensorial
% product of two 1D meshes. It is used when a non-obtuse triangular
% mesh is needed and in the self-consistent Monte Carlo device
% simulator. The following sheme shows the local numeration of
% nodes and edges:
%
%          in1 ---ie2--- in2
%             |         |
%             |         |
%            ie3       ie1
%             |         |
%             |         |
%          in4 ---ie4--- in3
%
%
% mesh.nn            number of nodes
% mesh.nt            number of triangles
% mesh.nr            number of rectangles
% mesh.ne            number of edges
% mesh.nrx           number of divisions along x
% mesh.nry           number of divisions along y
% mesh.nnx           number of points along x
% mesh.nny           number of points along y
% mesh.node:         Point Matrix (2,nn). The first and second rows contain
%                    x- and y-coordinates of the points in the mesh
% mesh.triangle:     Triangle Matrix (7,nt). Rows 1,2,3 contain indices to the corner points
%                    (nodes at vertices), given in counter clockwise order,
%                    and row 4 contains the subdomain numbers. Rows 5,6,7 contain indices to
%                    the three edges of the triangle.
% mesh.rect:         Rectangle Matrix (13,nr). Rows 1,2,3,4 contain indices to the corner points
%                    (nodes at vertices), given in clockwise order (see scheme above),
%                    and row 5 contains the subdomain numbers. Rows 6,7,8,9 contain
%                    indices to the edges of the rectangle. Rows 10,11,12,13
%                    contain indices to the neighbouring rectangles (used in FBMC).
% mesh.edge:         in the edge matrix E, the first and second rows contain
%                    indices of the starting and ending point,
%                    the third row contains the boundary segment, and
%                    the fourth row contains the boundary condition associated
%                    to the boundary segment (Monte Carlo):
%                       0 --> the particle is transmitted (black)
%                       1 --> semiconductor to dielectric, specular reflection (red)
%                       2 --> semiconductor to contact, absorption (blu)
%                       3 --> semiconductor to dielectric, periodic boundary condition
% mesh.contact:      (1,nn) for each node belonging to ith-contact
%                    stores index i, 0 otherwise
% mesh.edgemat       connectivity matrix
% mesh.dop:          net doping (1,nn), 1/cm^3
% mesh.xgrid         1D grid along x-direction, m
% mesh.ygrid         1D grid along y-direction, m
% mesh.lmean         average triangle side length, m
% mesh.inpec         mark points on perfect electric conductors (PEC)
% mesh.iepec         mark edges on perfect electric conductors (PEC)
% mesh.inpmc         mark points on perfect magnetic conductors (PMC)
% mesh.iepmc         mark edges on perfect magnetic conductors (PMC)
% mesh.ethj(1,:)     points lying on heterojunctions (side 1)
% mesh.ethj(2,:)     points lying on heterojunctions (side 2)
% mesh.ethj(3,:)     length of edges lying on heterojunctions, m
% mesh.point         triangle mid-points coordinates, m
% mesh.midpoint      rectangle mid-points coordinates, m
% mesh.contact_nadj  number of cells adajcent to contacts
% mesh.contact_seg   metal-semiconductor interfaces, used for Schottky contacts
% mesh.iq            mark points belonging to semiconductor regions
% mesh.gradx         sparse matrix (nn,nt), gradient along x
% mesh.grady         sparse matrix (nn,nt), gradient along y
% mesh.Se            array (nt), area associated to each triangle
% mesh.Ae            array (nn), area associated to each node
%
%
%
%
%
%
% Input variables:
%
% geom.dg      GEOMETRY DESCRIPTION MATRIX
%
%       The  Geometry Description matrix GD describes solid model that
%       you draw in the PDETOOL GUI.
%
%       Each column in the Geometry Description matrix corresponds to an
%       object in the solid geometry model. Four types of solid objects are
%       supported.  The object type is specified in row one:
%
%       1. For the circle solid, row one contains 1, the second and
%       third row contain the center x- and y-coordinates respectively.
%       Row four contains the radius of the circle.
%       2. For a polygon solid, row one contains 2, and the second row
%       contains the number, N, of line segments in the boundary.
%       The following N rows contain the x-coordinates of the starting points
%       of the edges, and the following N rows contain the y-coordinates of the
%       starting points of the edges.
%       3. For a rectangle solid, row one contains 3. The format is otherwise
%       identical to the polygon format.
%       4. For an ellipse solid, row one contains 4, the second and
%       third row contain the center x- and y-coordinates respectively. Row
%       four and five contain the major and minor axes of the ellipse.
%       The rotational angle of the ellipse is stored in row six.
%
% geom.dl      DECOMPOSED GEOMETRY MATRIX (geom.dl=decsg(geom.gd))
%
%       The Decomposed Geometry matrix DL contains a representation of
%       the minimal regions that have been constructed by the DECSG algorithm.
%       Each edge segment of the minimal regions correspond to a column
%       in DL. In each such column rows two and three contain the starting
%       and ending x-coordinate, and rows four and five the corresponding
%       y-coordinate. Rows six and seven contain left and right minimal region
%       numbers with respect to the direction induced by the start and end
%       points (counter clockwise direction on circle and ellipse segments).
%       There are three types of possible edge segments in a minimal regions:
%
%       1. For circle edge segments row one is 1. Rows eight and nine
%       contain the coordinates of the center of the circle. Row 10 contains
%       the radius.
%       2. For line edge segments row one is 2.
%       3. For ellipse edge segments row one is 4. Rows eight and nine
%       contain the coordinates of the center of the ellipse. Rows 10
%       and 11 contain the major and minor axes of the ellipse
%       respectively. The rotational angle of the ellipse is stored in
%       row 12.
%
% geom.gdm     DECOMPOSED GEOMETRY MATRIX (geom.dgm=computedgm(geom))
%
%===============================================================================================100
% command preprocessing

%'entro', keyboard
if(not(isfield(geom,'mobility_dop'))), geom.mobility_dop='none'; end
if(not(isfield(geom,'ethjflg'))), geom.ethjflg=0; end
if(not(isfield(geom,'dop_t'))), geom.dop_t=[]; end
%===============================================================================================100
I=find(not(cellfun('isempty',geom.X))); J=find(not(cellfun('isempty',geom.Y)));
xb=evalgeom(geom,geom.X(I)); yb=evalgeom(geom,geom.Y(J));
Nx=sum(geom.div_x)+1; Ny=sum(geom.div_y)+1; nrx=Nx-1; nry=Ny-1; nr=nrx*nry; nn=Nx*Ny;
x=zeros(1,Nx); y=zeros(1,Ny);
index=1;
for k=1:length(geom.div_x),
    x(index:sum(geom.div_x(1:k))+1)=linspace(xb(k),xb(k+1),geom.div_x(k)+1);
    index=index+geom.div_x(k); end
index=1;
for k=1:length(geom.div_y),
    y(index:sum(geom.div_y(1:k))+1)=linspace(yb(k),yb(k+1),geom.div_y(k)+1);
    index=index+geom.div_y(k); end
mesh.xgrid=x; mesh.ygrid=y;
%===============================================================================================100
% compute node matrix
[x,y]=meshgrid(x,y);
ix=1:Nx; iy=1:Ny; [ix,iy]=meshgrid(ix,iy);
in=(ix-1)*Ny+iy; % nodal matrix index
in1=reshape(in(2:Ny,1:Nx-1)',1,nr);
in2=reshape(in(2:Ny,2:Nx)',1,nr);
in3=reshape(in(1:Ny-1,2:Nx)',1,nr);
in4=reshape(in(1:Ny-1,1:Nx-1)',1,nr);
x1=x(in1); x2=x(in2); x3=x(in3); x4=x(in4);
y1=y(in1); y2=y(in2); y3=y(in3); y4=y(in4);
l1=y2-y3; l2=x2-x1; l3=y1-y4; l4=x3-x4;
mesh.lmean=sum(l1+l2+l3+l4)/(4*nr);
xc=(x1+x2)./2; yc=(y2+y3)./2;
x=reshape(x,1,nn); y=reshape(y,1,nn);
mesh.node=[x; y];
%===============================================================================================100
% compute rect matrix
region=zeros(1,nr);
nreg=max([geom.dgm(6,:) geom.dgm(7,:)]);
for sdl=1:nreg;
    [xv,yv]=extractpatch(geom.dgm,sdl);
    I=inpolygon(xc,yc,xv,yv); I=find(I); region(I)=sdl; end
iemat=reshape(cumsum(ones(1,(Ny-1)*Nx)),Nx,Ny-1)';
ie1=reshape(iemat(:,2:Nx)',1,nr);
ie3=reshape(iemat(:,1:Nx-1)',1,nr);
iemat=reshape(cumsum(ones(1,Ny*(Nx-1))),Nx-1,Ny)';
iemat=iemat+(Ny-1)*Nx;
ie4=reshape(iemat(1:Ny-1,:)',1,nr);
ie2=reshape(iemat(2:Ny,:)',1,nr);
irmat=zeros(Ny+1,Nx+1);
irmat(2:Ny,2:Nx)=reshape(cumsum(ones(1,nr)),Nx-1,Ny-1)';
ir1=reshape(irmat(2:Ny,3:Nx+1)',1,nr);
ir3=reshape(irmat(2:Ny,1:Nx-1)',1,nr);
ir4=reshape(irmat(1:Ny-1,2:Nx)',1,nr);
ir2=reshape(irmat(3:Ny+1,2:Nx)',1,nr);
mesh.rect=[in1; in2; in3; in4; region; ie1; ie2; ie3; ie4; ir1; ir2; ir3; ir4];
mesh.nr=nr; mesh.nrx=nrx; mesh.nry=nry;
mesh.nn=nn; mesh.nnx=Nx; mesh.nny=Ny;
mesh.ne=(Nx-1)*Ny+Nx*(Ny-1);
%===============================================================================================100
% compute edge matrix
mesh=f_edge(geom,mesh);
%===============================================================================================100
% duplicate points at heterojunctions
if(geom.ethjflg), mesh=f_ethj(geom,mesh); nn=mesh.nn; end;
% kludge: f_ethj should also update the edge matrix!
%===============================================================================================100
% compute contact data
mesh=f_contact(geom,mesh);
%===============================================================================================100
% drift mesh (used to track particles in Monte Carlo)
[mesh2]=rectmesh2(geom);
mesh.node2=mesh2.node; mesh.rect2=mesh2.rect; mesh.edge2=mesh2.edge;
mesh.nr2=mesh2.nr; mesh.nrx2=mesh2.nrx; mesh.nry2=mesh2.nry;
mesh.nn2=mesh2.nn; mesh.nnx2=mesh2.nnx; mesh.nny2=mesh2.nny;
mesh.ne2=mesh2.ne;
mesh.xgrid2=mesh2.xgrid;
mesh.ygrid2=mesh2.ygrid;
%===============================================================================================100
% compute triangle matrix
mesh.triangle = [mesh.rect(1,:) mesh.rect(3,:); ...
    mesh.rect(4,:) mesh.rect(2,:); ...
    mesh.rect(2,:) mesh.rect(4,:); ...
    mesh.rect(5,:) mesh.rect(5,:)];
nt = 2*mesh.nr; mesh.nt=nt;
in1=mesh.triangle(1,:); in2=mesh.triangle(2,:); in3=mesh.triangle(3,:);
% make a connectivity matrix (if nodes i and j are connected A(i,j)>0 , otherwise A(i,j)=0)
Amat=sparse([in1 in2 in3],[in2 in3 in1],1,nn,nn); Amat=Amat+Amat';
[ieb,iee]=find(Amat); I=ieb<iee; ieb=ieb(I); iee=iee(I);
nedges=length(ieb); mesh.edgemat=sparse(ieb,iee,1:nedges,nn,nn);
mesh.edgemat=mesh.edgemat-mesh.edgemat';
for t=1:nt;
    mesh.triangle(5,t)=mesh.edgemat(in2(t),in3(t));
    mesh.triangle(6,t)=mesh.edgemat(in3(t),in1(t));
    mesh.triangle(7,t)=mesh.edgemat(in1(t),in2(t)); end
mesh.nedges=nedges;
%===============================================================================================100
% find nodes on perfect electric walls
mesh.inpec=logical(zeros(1,mesh.nn));
for ii=1:length(geom.bspec);
    mesh.inpec(mesh.edge(1,mesh.edge(3,:)==geom.bspec(ii)))=1;
    mesh.inpec(mesh.edge(2,mesh.edge(3,:)==geom.bspec(ii)))=1; end
% find nodes on perfect magnetic walls
mesh.inpmc=logical(zeros(1,mesh.nn));
for ii=1:length(geom.bspmc);
    mesh.inpmc(mesh.edge(1,mesh.edge(3,:)==geom.bspmc(ii)))=1;
    mesh.inpmc(mesh.edge(1,mesh.edge(3,:)==geom.bspmc(ii)))=1; end
%===============================================================================================100
if(not(geom.ethjflg)),
    % find edges on perfect electric walls
    mesh.iepec=logical(zeros(1,mesh.nedges));
    inpec=mesh.edge(1:2,ismember(mesh.edge(3,:),geom.bspec)); iepec=[];
    for ii=1:size(inpec,2); iepec(ii)=abs(mesh.edgemat(inpec(1,ii),inpec(2,ii))); end
    mesh.iepec(iepec)=1;
    % find edges on perfect magnetic walls
    mesh.iepmc=logical(zeros(1,mesh.nedges));
    inpmc=mesh.edge(1:2,ismember(mesh.edge(3,:),geom.bspmc)); iepmc=[];
    for ii=1:size(inpmc,2); iepmc(ii)=abs(mesh.edgemat(inpmc(1,ii),inpmc(2,ii))); end
    mesh.iepmc(iepmc)=1;
end
%===============================================================================================100
% define doping profile
mesh=f_dop(geom,mesh);
if(not(isempty(geom.function_dop))),
    eval(['mesh=',geom.function_dop,'(geom,mesh);']), end
mesh.dop_d_t=pdeintrp(mesh.node,mesh.triangle(1:4,:),mesh.dop_d.');
mesh.dop_a_t=pdeintrp(mesh.node,mesh.triangle(1:4,:),mesh.dop_a.');
%===============================================================================================100
% define trap concentrations
if(not(isempty(geom.dop_t))), [mesh]=f_trap(geom,mesh); end
%===============================================================================================100
% compute the midpoints for the triangular mesh
in1=mesh.triangle(1,:); in2=mesh.triangle(2,:); in3=mesh.triangle(3,:);
x1=mesh.node(1,in1); x2=mesh.node(1,in2); x3=mesh.node(1,in3);
y1=mesh.node(2,in1); y2=mesh.node(2,in2); y3=mesh.node(2,in3);
mesh.point=[(x1+x2+x3)/3; (y1+y2+y3)/3];
% mark the points belonging to semiconductor regions
it=ismember(mesh.triangle(4,:),find(geom.semiconductor)); % triangles in semiconductor regions
iq=false(1,nn); iq([in1(it) in2(it) in3(it)])=1; mesh.iq=iq;
% compute gradient operators
r3x=x2-x1; r1x=x3-x2; r2x=x1-x3;
r3y=-y2+y1; r1y=-y3+y2; r2y=-y1+y3;
f11xx=r1x.*r1x; f12xx=r1x.*r2x; f13xx=r1x.*r3x;
f22xx=r2x.*r2x; f23xx=r2x.*r3x; f31xx=f13xx; f33xx=r3x.*r3x;
f11yy=r1y.*r1y; f12yy=r1y.*r2y; f13yy=r1y.*r3y;
f22yy=r2y.*r2y; f23yy=r2y.*r3y; f31yy=f13yy; f33yy=r3y.*r3y;
f23xy=r2x.*r3y; f32xy=r3x.*r2y;
% triangle side lengths
xm1=1/2*(x2+x3); xm2=1/2*(x3+x1); xm3=1/2*(x1+x2);
ym1=1/2*(y2+y3); ym2=1/2*(y3+y1); ym3=1/2*(y1+y2);
% split points
xcc=1./2.*(x2.*f32xy+x1.*f32xy-y2.*f23yy-x3.*f23xy-x1.*f23xy+y3.*f23yy)./(f32xy-f23xy);
ycc=1./2.*(x2.*f23xx-y2.*f23xy-y1.*f23xy-x3.*f23xx+y3.*f32xy+y1.*f32xy)./(f32xy-f23xy);
% triangle mid-points
s1=sqrt((xm1-xcc).^2+(ym1-ycc).^2)*0;
s2=sqrt((xm2-xcc).^2+(ym2-ycc).^2);
s3=sqrt((xm3-xcc).^2+(ym3-ycc).^2);
Se=abs(r2x.*r1y-r2y.*r1x)./2; % areas associated to triangles
mesh.gradx=1/2*sparse([in1 in2 in3],[1:nt 1:nt 1:nt],[r1y./Se r2y./Se r3y./Se],nn,nt);
mesh.grady=1/2*sparse([in1 in2 in3],[1:nt 1:nt 1:nt],[r1x./Se r2x./Se r3x./Se],nn,nt);
% set triangle area to zero outside semiconductors
it=not(ismember(mesh.triangle(4,:),find(geom.semiconductor)));
Se(it)=0; mesh.Se=Se;
% compute areas associated to nodes in semiconductor regions
Se1=1/2*(xm3.*ycc-xcc.*ym3+x1.*ym3-y1.*xm3+xcc.*ym2-xm2.*ycc-x1.*ym2+y1.*xm2);
Se2=1/2*(xm1.*ycc-xcc.*ym1+x2.*ym1-y2.*xm1+xcc.*ym3-xm3.*ycc-x2.*ym3+y2.*xm3);
Se3=1/2*(xm2.*ycc-xcc.*ym2+x3.*ym2-y3.*xm2+xcc.*ym1-xm1.*ycc-x3.*ym1+y3.*xm1);
% set triangle area to zero outside semiconductors
it=not(ismember(mesh.triangle(4,:),find(geom.semiconductor)));
Se1(it)=0; Se2(it)=0; Se3(it)=0;
mesh.Ae=full(sparse(1,[in1 in2 in3],[Se1 Se2 Se3],1,nn));
%==============================================================================================100
% compute the midpoints for the rectangular mesh (Monte Carlo)
in1=mesh.rect(1,:); in2=mesh.rect(2,:);
in3=mesh.rect(3,:); in4=mesh.rect(4,:);
x1=mesh.node(1,in1); x2=mesh.node(1,in2);
y1=mesh.node(2,in1); y4=mesh.node(2,in4);
mesh.midpoint=[(x1+x2)/2; (y1+y4)/2];
%
%==============================================================================================100
% tagging quantum well (QW) cells
%==============================================================================================100
if(isfield(geom,'QWorientation'))
    % initialization of mesh in matrix form
    xx=mesh.node(1,:); yy=mesh.node(2,:);
    mesh.X=reshape(xx,mesh.nny,mesh.nnx); mesh.Y=reshape(yy,mesh.nny,mesh.nnx);
    mesh.IN=reshape(1:mesh.nn,mesh.nny,mesh.nnx);
    %
    % initialization of variables to be saved
    NQW=0;
    %
    for indCell=1:length(geom.label)
        VerQW=strfind(geom.label{indCell},'qw'); % looking for QW
        VerQWPassiv=strfind(geom.label{indCell},'PassivationQW'); % looking for passivation out of the QW
        if(not(isempty(VerQW)))
            NQW=NQW+1;
            % index of the quantum well from Debernardi label
            indMQW=str2num(geom.label{indCell}(3));
            %
            % common quantities for QWs
            ITrQW=find(mesh.triangle(4,:)==indCell);
            IQW=unique([mesh.triangle(1,ITrQW),mesh.triangle(2,ITrQW),mesh.triangle(3,ITrQW)]);
            %
            % mesh.IQW=[mesh.IQW,IQW]; % IQW is for 3D Ccap terms
            % for x (i.e. also rho) oriented QWs
            WCav=max(mesh.node(2,IQW))-min(mesh.node(2,IQW));
            %
            % computing mesh node closest to center of the QW
            % it is assumed to have x(or rho)-oriented quantum wells
            yQWMean=mean(mesh.node(2,IQW));
            [~,iQW]=min(abs(yQWMean-mesh.ygrid));
            %
            mesh.yMQW{indMQW}=mesh.node(2,iQW);
            inMQW = find(ismember(mesh.node(2,:),mesh.yMQW{indMQW}));
            min_xQW = min(mesh.node(1,IQW));
            max_xQW = max(mesh.node(1,IQW));
            inMQW = inMQW(mesh.node(1,inMQW) >= min_xQW & mesh.node(1,inMQW) <= max_xQW);
            %
            mesh.vWMQW{indMQW}=WCav;
            mesh.MQWcell{indMQW}=indCell;
            %
            if(not(isfield(mesh,'inMQW')))
                % Excluding passivation
                mesh.inMQW{indMQW}=inMQW;
                mesh.IMQW{indMQW}=IQW;
                mesh.ITrMQW{indMQW}=ITrQW;
                % Including passivation
                mesh.inMQWP{indMQW}=inMQW;
                mesh.IMQWP{indMQW}=IQW;
                mesh.ITrMQWP{indMQW}=ITrQW;
            else
                % Excluding passivation
                mesh.inMQW{indMQW}=unique([mesh.inMQW{indMQW},inMQW]);
                mesh.IMQW{indMQW}=unique([mesh.IMQW{indMQW},IQW]);
                mesh.ITrMQW{indMQW}=[mesh.ITrMQW{indMQW},ITrQW];
                % Including passivation
                mesh.inMQWP{indMQW}=unique([mesh.inMQW{indMQW},inMQW]);
                mesh.IMQWP{indMQW}=[mesh.IMQW{indMQW},IQW];
                mesh.ITrMQWP{indMQW}=[mesh.ITrMQW{indMQW},ITrQW];
            end
            %
        end
        if(not(isempty(VerQWPassiv)))
            % index of the quantum well from Debernardi label
            indMQW=str2num(geom.label{indCell}(14));
            %
            % common quantities for QWs
            ITrQW=find(mesh.triangle(4,:)==indCell);
            IQW=unique([mesh.triangle(1,ITrQW),mesh.triangle(2,ITrQW),mesh.triangle(3,ITrQW)]);
            %
            % mesh.IQW=[mesh.IQW,IQW]; % IQW is for 3D Ccap terms
            % for x (i.e. also rho) oriented QWs
            WCav=max(mesh.node(2,IQW))-min(mesh.node(2,IQW));
            %
            % computing mesh node closest to center of the QW
            % it is assumed to have x(or rho)-oriented quantum wells
            yQWMean=mean(mesh.node(2,IQW));
            [~,iQW]=min(abs(yQWMean-mesh.ygrid));
            %
            inMQW = find(ismember(mesh.node(2,:),mesh.yMQW{indMQW}));
            min_xQW = min(mesh.node(1,IQW));
            max_xQW = max(mesh.node(1,IQW));
            inMQW = inMQW(mesh.node(1,inMQW) >= min_xQW & mesh.node(1,inMQW) <= max_xQW);
            %
            if(not(isfield(mesh,'inMQWP')))
                mesh.inMQWP{indMQW}=inMQW;
                mesh.IMQWP{indMQW}=IQW;
                mesh.ITrMQWP{indMQW}=ITrQW;
            else
                mesh.inMQWP{indMQW}=unique([mesh.inMQW{indMQW},inMQW]);
                mesh.IMQWP{indMQW}=unique([mesh.IMQW{indMQW},IQW]);
                mesh.ITrMQWP{indMQW}=[mesh.ITrMQW{indMQW},ITrQW];
            end
        end
    end
    %
    mesh.NMQW=length(mesh.vWMQW);
    %
    %
    % Loading the LUT from the k.p subband solver
    eval(['load ',geom.GLUTm])
    %
    for indQW=1:mesh.NMQW
        %
        mesh.nnxQW{indQW}=length(mesh.inMQW{indQW}); % number of QW x nodes
        % Saving k.p eigenvalues and eigenfunctions
        mesh.lambda_C{indQW}=Ban.SBC(:,1)-meshQW.Eg;
        mesh.lambda_V{indQW}=Ban.SBV(:,1);
        Psi_C2_QW=abs(Ban.XVC).^2.'; % conduction band
        Psi_V2_QW=abs(Ban.XVV).^2.'; % valence band
        %
        % Loading QW geometry, material parameters and mesh from k.p
        y=(meshQW.x-meshQW.L/2)*1e2; % shift mesh origin in QW center, cm
        WQW=meshQW.Lz*1e2; % quantum well width from k.p, cm
        NC=size(Ban.SBC,1);
        mesh.meffnMQW{indQW}=ones(NC,1)*meshQW.meffn;
        mesh.meffpMQW{indQW}=MefH.';
        mesh.xmol_barrierLUT=meshQW.xmol_barrier;
        mesh.xmol_wellLUT=meshQW.xmol_well;
        mesh.QcLUT=meshQW.Qc;
        mesh.DeltaEgLUT=meshQW.DeltaEg;
        mesh.DeltaEcQW{indQW}=meshQW.DeltaEg*meshQW.xmol_barrier*meshQW.Qc;
        mesh.DeltaEvQW{indQW}=meshQW.DeltaEg*meshQW.xmol_barrier*(1-meshQW.Qc);
        %
        % Identifying regions before, in and after the quantum well
        indQWLeft=find(y<-WQW/2); % points before the quantum well
        % indQWCenter=find(y>=-WQW/2 & y<=WQW/2); % quantum well points
        indQWRight=find(y>WQW/2); % points after the quantum well
        %
        % Indexes for interpolation of the eigenfunctions out of the QW
        indInterpLeft=indQWLeft(end-10:end-3);
        indInterpRight=indQWRight(3:10);
        %
        % Performing a linear interpolation of the logarithm of each
        % conduction-like eigenfunction, to describe the exponential decay
        % out of the quantum well
        for indC=1:size(Psi_C2_QW,1)
            CoeffLeftC(indC,:)=polyfit(y(indInterpLeft),log(abs(Psi_C2_QW(indC,indInterpLeft))),1);
            CoeffRightC(indC,:)=polyfit(y(indInterpRight),log(abs(Psi_C2_QW(indC,indInterpRight))),1);
        end
        %
        % Performing a linear interpolation of the logarithm of each
        % valence-like eigenfunction, to describe the exponential decay
        % out of the quantum well
        for indV=1:size(Psi_V2_QW,1)
            CoeffLeftV(indV,:)=polyfit(y(indInterpLeft),log(abs(Psi_V2_QW(indV,indInterpLeft))),1);
            CoeffRightV(indV,:)=polyfit(y(indInterpRight),log(abs(Psi_V2_QW(indV,indInterpRight))),1);
        end
        %
        % Identifying regions before, in and after the quantum well, in
        % drift-diffusion mesh
        yQW=mesh.yMQW{indQW};
        indDDLeft=find(mesh.ygrid<(yQW-WQW/2));
        indDDCenter=find(mesh.ygrid>=(yQW-WQW/2) & mesh.ygrid<=(yQW+WQW/2));
        indDDRight=find(mesh.ygrid>(yQW+WQW/2));
        yLeft=mesh.ygrid(indDDLeft);
        yCenter=mesh.ygrid(indDDCenter);
        yRight=mesh.ygrid(indDDRight);
        %
        % Interpolation of conduction eigenfunctions in DD nodes
        for indC=1:size(Psi_C2_QW,1)
            Psi_C2_DD(indC,indDDLeft)=exp(polyval(CoeffLeftC(indC,:),yLeft-yQW));
            Psi_C2_DD(indC,indDDRight)=exp(polyval(CoeffRightC(indC,:),yRight-yQW));
            Psi_C2_DD(indC,indDDCenter)=interp1(y,Psi_C2_QW(indC,:),yCenter-yQW);
            NormC2=trapz(mesh.ygrid,Psi_C2_DD(indC,:));
            Psi_C2_DD(indC,:)=Psi_C2_DD(indC,:)/NormC2;
        end
        %
        % Interpolation of valence eigenfunctions in DD nodes
        for indV=1:size(Psi_V2_QW,1)
            Psi_V2_DD(indV,indDDLeft)=exp(polyval(CoeffLeftV(indV,:),yLeft-yQW));
            Psi_V2_DD(indV,indDDRight)=exp(polyval(CoeffRightV(indV,:),yRight-yQW));
            Psi_V2_DD(indV,indDDCenter)=interp1(y,Psi_V2_QW(indV,:),yCenter-yQW);
            NormV2=trapz(mesh.ygrid,Psi_V2_DD(indV,:));
            Psi_V2_DD(indV,:)=Psi_V2_DD(indV,:)/NormV2;
        end
        %
        mesh.Psi_C2{indQW}=Psi_C2_DD;
        mesh.Psi_V2{indQW}=Psi_V2_DD;
        
        % % DEBUG PLOTS
%         figure,hold on,plot(y,Psi_C2_QW/1e2,'-'),set(gca,'yscale','log')
%         ax = gca;
%         ax.ColorOrderIndex = 1;
%         plot(mesh.ygrid-yQW,Psi_C2_DD,'o'),set(gca,'yscale','log')
%         
%         figure,hold on,plot(y,Psi_V2_QW/1e2,'-'),set(gca,'yscale','log')
%         ax = gca;
%         ax.ColorOrderIndex = 1;
%         plot(mesh.ygrid-yQW,Psi_V2_DD,'o'),set(gca,'yscale','log')
        
        % mesh.lambda_C{1}
        % mesh.lambda_V{1}
    end
    % Finding cavity indexes for bulk Auger expulsion
    flagCav=-1; % -1 indicates NOT to save indices
    flagCavOld=flagCav; % old cavity flag
    for indCell=1:length(geom.label)
        VerCav=strfind(geom.label{indCell},'Cav'); % looking for QW
        % save current (that becomes old) flagCav value
        flagCavOld = flagCav;
        if(not(isempty(VerCav)) & VerCav==1)
            % toggle flagCav, so that the code starts looking for cavity
            % indexes
            flagCav = - flagCav; 
        end
        % The following operations should be performed for all the cells
        % within two "Cav" labels, including the two Cavs. For this reason,
        % it is useful to save both flagCav and flagCavOld. The two
        % conditions for entering the following instruction set are either
        % - having flagCav==1, i.e., meeting the first 'Cav' label
        % or
        % - having flagCav==-1, and flagCavOld==+1, i.e., having met the
        % last 'Cav' label
        if((flagCav==1 | (flagCav==-1 & flagCavOld==1)) & strcmp(geom.material(indCell),'Polyamide')==0)
            ITrCav=find(mesh.triangle(4,:)==indCell);
            ICAV=unique([mesh.triangle(1,ITrCav),mesh.triangle(2,ITrCav),mesh.triangle(3,ITrCav)]);
            %
            %
            if(not(isfield(mesh,'ICAV')))
                mesh.CAVcell=indCell;
                mesh.ICAV=ICAV;
                mesh.ITrCAV=ITrCav;
            else
                mesh.ICAV=unique([mesh.ICAV,ICAV]);
                mesh.ITrCAV=unique([mesh.ITrCAV,ITrCav]);
                if(isempty(mesh.CAVcell))
                    mesh.CAVcell=indCell;
                    mesh.ICAV=ICAV;
                    mesh.ITrCAV=ITrCav;
                end
            end
        end
    end
    % for x (i.e. also rho) oriented QWs
    WCav=max(mesh.node(2,mesh.ICAV))-min(mesh.node(2,mesh.ICAV));
    % one-time defined quantities (put here! not out!)
    mesh.WCav=WCav;
else
    mesh.NMQW=0;
end % if(isfield(mode,'QWorientation'))

if(isfield(geom,'iNEGF'))
    %
%     for indCell=1:(length(geom.label)/length(geom.div_x))
    for indCell=1:length(geom.label)
        VerBTJ=strfind(geom.label{indCell},'TJ'); % looking for BTJ layers
        VerBTJPassiv=strfind(geom.label{indCell},'PassivationTunn'); % looking for passivation out of the QW
        %
        if(not(isempty(VerBTJ)))% && not(strcmp(geom.material{indCell},'AlOx')))
            %
            indBTJ_layer=str2double(geom.label{indCell}(end));
            indBTJ=str2double(geom.label{indCell}(end-2));
            %
            % common quantities for BTJs
            ITrBTJ=find(mesh.triangle(4,:)==indCell);
%             IBTJ=ITrBTJ;
            IBTJ=unique([mesh.triangle(1,ITrBTJ),mesh.triangle(2,ITrBTJ),mesh.triangle(3,ITrBTJ)]);
            %
            % computing mesh node closest to center of the QW
            % it is assumed to have x(or rho)-oriented quantum wells
            WBTJ=max(mesh.node(2,IBTJ))-min(mesh.node(2,IBTJ));
            %
            yBTJLeft=min(mesh.node(2,IBTJ));
            yBTJRight=max(mesh.node(2,IBTJ));
            yBTJMean=mean(mesh.node(2,IBTJ));
            [~,iBTJ]=min(abs(yBTJMean-mesh.ygrid));
            [~,iLeftBTJ]=min(abs(yBTJLeft-mesh.ygrid));
            [~,iRightBTJ]=min(abs(yBTJRight-mesh.ygrid));
            %
            mesh.yBTJ{indBTJ,indBTJ_layer}=mesh.node(2,iBTJ);
            mesh.yLeftBTJ{indBTJ,indBTJ_layer}=mesh.node(2,iLeftBTJ);
            mesh.yRightBTJ{indBTJ,indBTJ_layer}=mesh.node(2,iRightBTJ);
            mesh.iLeftBTJ{indBTJ,indBTJ_layer}=iLeftBTJ;
            mesh.iRightBTJ{indBTJ,indBTJ_layer}=iRightBTJ;
            %
            inBTJ = find(ismember(mesh.node(2,:),mesh.yBTJ{indBTJ,indBTJ_layer}));
            min_xBTJ = min(mesh.node(1,IBTJ));
            max_xBTJ = max(mesh.node(1,IBTJ));
            inBTJ = inBTJ(mesh.node(1,inBTJ) >= min_xBTJ & mesh.node(1,inBTJ) <= max_xBTJ);
            %
            iLBTJ = find(ismember(mesh.node(2,:),mesh.yLeftBTJ{indBTJ,indBTJ_layer}));
            iLBTJ = iLBTJ(mesh.node(1,iLBTJ) >= min_xBTJ & mesh.node(1,iLBTJ) <= max_xBTJ);
            %
            iRBTJ = find(ismember(mesh.node(2,:),mesh.yRightBTJ{indBTJ,indBTJ_layer}));
            iRBTJ = iRBTJ(mesh.node(1,iRBTJ) >= min_xBTJ & mesh.node(1,iRBTJ) <= max_xBTJ);
            %
            mesh.vBTJ{indBTJ,indBTJ_layer}=WBTJ;
            mesh.BTJcell{indBTJ,indBTJ_layer}=indCell;
            %
            if(not(isfield(mesh,'inBTJ')))
                % Excluding passivation
                mesh.inBTJ{indBTJ,indBTJ_layer}=inBTJ;
                mesh.iLBTJ{indBTJ,indBTJ_layer}=iLBTJ;
                mesh.iRBTJ{indBTJ,indBTJ_layer}=iRBTJ;
                mesh.IBTJ{indBTJ,indBTJ_layer}=IBTJ;
                mesh.ITrBTJ{indBTJ,indBTJ_layer}=ITrBTJ;
                % Including passivation
                mesh.inBTJP{indBTJ,indBTJ_layer}=inBTJ;
                mesh.iLBTJP{indBTJ,indBTJ_layer}=iLBTJ;
                mesh.iRBTJP{indBTJ,indBTJ_layer}=iRBTJ;
                mesh.IBTJP{indBTJ,indBTJ_layer}=IBTJ;
                mesh.ITrBTJP{indBTJ,indBTJ_layer}=ITrBTJ;
            else
                % Excluding passivation
                mesh.inBTJ{indBTJ,indBTJ_layer}=unique([mesh.inBTJ{indBTJ,indBTJ_layer},inBTJ]);
                mesh.iLBTJ{indBTJ,indBTJ_layer}=unique([mesh.iLBTJ{indBTJ,indBTJ_layer},iLBTJ]);
                mesh.iRBTJ{indBTJ,indBTJ_layer}=unique([mesh.iRBTJ{indBTJ,indBTJ_layer},iRBTJ]);
                mesh.IBTJ{indBTJ,indBTJ_layer}=unique([mesh.IBTJ{indBTJ,indBTJ_layer},IBTJ]);
                mesh.ITrBTJ{indBTJ,indBTJ_layer}=[mesh.ITrBTJ{indBTJ,indBTJ_layer},ITrBTJ];
                % Including passivation
                mesh.inBTJP{indBTJ,indBTJ_layer}=unique([mesh.inBTJ{indBTJ,indBTJ_layer},inBTJ]);
                mesh.iLBTJP{indBTJ,indBTJ_layer}=unique([mesh.iLBTJ{indBTJ,indBTJ_layer},iLBTJ]);
                mesh.iRBTJP{indBTJ,indBTJ_layer}=unique([mesh.iRBTJ{indBTJ,indBTJ_layer},iRBTJ]);
                mesh.IBTJP{indBTJ,indBTJ_layer}=[mesh.IBTJ{indBTJ,indBTJ_layer},IBTJ];
                mesh.ITrBTJP{indBTJ,indBTJ_layer}=[mesh.ITrBTJ{indBTJ,indBTJ_layer},ITrBTJ];
            end
        end
        if(not(isempty(VerBTJPassiv)))
            % index of the quantum well from Debernardi label
            indBTJ_layer=str2double(geom.label{indCell}(end));
            %
            % common quantities for QWs
            ITrBTJ=find(mesh.triangle(4,:)==indCell);
%             IBTJ=ITrBTJ;
            IBTJ=unique([mesh.triangle(1,ITrBTJ),mesh.triangle(2,ITrBTJ),mesh.triangle(3,ITrBTJ)]);
            %
            % mesh.IQW=[mesh.IQW,IQW]; % IQW is for 3D Ccap terms
            % for x (i.e. also rho) oriented QWs
            WBTJ=max(mesh.node(2,IBTJ))-min(mesh.node(2,IBTJ));
            %
            % computing mesh node closest to center of the QW
            % it is assumed to have x(or rho)-oriented quantum wells
            yBTJMean=mean(mesh.node(2,IBTJ));
            [~,iBTJ]=min(abs(yBTJMean-mesh.ygrid));
            %
            inBTJ = find(ismember(mesh.node(2,:),mesh.yBTJ{indBTJ,indBTJ_layer}));
            min_xBTJ = min(mesh.node(1,IBTJ));
            max_xBTJ = max(mesh.node(1,IBTJ));
            inBTJ = inBTJ(mesh.node(1,inBTJ) >= min_xBTJ & mesh.node(1,inBTJ) <= max_xBTJ);
            %
            iLBTJ = find(ismember(mesh.node(2,:),mesh.yLeftBTJ{indBTJ,indBTJ_layer}));
            iLBTJ = iLBTJ(mesh.node(1,iLBTJ) >= min_xBTJ & mesh.node(1,iLBTJ) <= max_xBTJ);
            %
            iRBTJ = find(ismember(mesh.node(2,:),mesh.yRightBTJ{indBTJ,indBTJ_layer}));
            iRBTJ = iRBTJ(mesh.node(1,iRBTJ) >= min_xBTJ & mesh.node(1,iRBTJ) <= max_xBTJ);
            %
            if(not(isfield(mesh,'inBTJP')))
                mesh.inBTJP{indBTJ,indBTJ_layer}=inBTJ;
                mesh.iLBTJP{indBTJ,indBTJ_layer}=iLBTJ;
                mesh.iRBTJP{indBTJ,indBTJ_layer}=iRBTJ;
                mesh.IBTJP{indBTJ,indBTJ_layer}=IBTJ;
                mesh.ITrBTJP{indBTJ,indBTJ_layer}=ITrBTJ;
            else
                mesh.inBTJP{indBTJ,indBTJ_layer}=unique([mesh.inBTJ{indBTJ,indBTJ_layer},inBTJ]);
                mesh.iLBTJP{indBTJ,indBTJ_layer}=unique([mesh.iLBTJ{indBTJ,indBTJ_layer},iLBTJ]);
                mesh.iRBTJP{indBTJ,indBTJ_layer}=unique([mesh.iRBTJ{indBTJ,indBTJ_layer},iRBTJ]);
                mesh.IBTJP{indBTJ,indBTJ_layer}=unique([mesh.IBTJ{indBTJ,indBTJ_layer},IBTJ]);
                mesh.ITrBTJP{indBTJ,indBTJ_layer}=[mesh.ITrBTJ{indBTJ,indBTJ_layer},ITrBTJ];
            end
        end        
    end
    %
    if isfield(mesh,'IBTJ')
        mesh.LBTJ=mesh.iLBTJ{end};  % left extreme points of the BTJ (needed for the derivatives)
        mesh.RBTJ=mesh.iRBTJ{1};    % right extreme points of the BTJ (needed for the derivatives)
    end
end
%
%***********************************************************************************************100
%
%
%***********************************************************************************************100
function [x_]=evalgeom(geom_,X_);
eval(char(geom_.par)')
for k_=1:length(X_); x_(k_)=eval(char(X_(k_))); end
%***********************************************************************************************100
%
%
%***********************************************************************************************100
function [xv,yv]=extractpatch(dl,sdl);
Ibs=ismember(dl(6,:),sdl)|ismember(dl(7,:),sdl);
x=dl(2:3,Ibs); y=dl(4:5,Ibs);
xv=x(:,1); yv=y(:,1); x(:,1)=NaN; y(:,1)=NaN;
for n=2:(size(x,2)-1)
    [dmin1,I1]=min((x(1,:)-xv(length(xv))).^2+(y(1,:)-yv(length(yv))).^2);
    [dmin2,I2]=min((x(2,:)-xv(length(xv))).^2+(y(2,:)-yv(length(yv))).^2);
    if(dmin1<dmin2), xv=[xv; x(2,I1)]; yv=[yv; y(2,I1)]; x(:,I1)=NaN; y(:,I1)=NaN;
    else, xv=[xv; x(1,I2)]; yv=[yv; y(1,I2)]; x(:,I2)=NaN; y(:,I2)=NaN; end
end
%***********************************************************************************************100
%
%
%**********************************************************************************************100
function [mesh]=f_contact(geom,mesh);
% contact matrix
mesh.contact=zeros(1,mesh.nn);
if(geom.nc==0), return, end
for sd=1:geom.nd;
    if(geom.reg2contact(sd)>0),
        ic=geom.reg2contact(sd);
        ii=ismember(mesh.rect(5,:),sd);
        mesh.contact(mesh.rect(1,ii))=ic;
        mesh.contact(mesh.rect(2,ii))=ic;
        mesh.contact(mesh.rect(3,ii))=ic;
        mesh.contact(mesh.rect(4,ii))=ic;
    end, end
%==============================================================================================100
% here we find the rectangles adjacent to contacts, used in FBMC
for ic=1:geom.nc;
    segment=find(ismember(geom.dgm(6,:),find(geom.reg2contact==ic)) | ...
        ismember(geom.dgm(7,:),find(geom.reg2contact==ic)));
    ie=ismember(mesh.edge(3,:),segment);
    in=zeros(1,mesh.nn); in(mesh.edge(1,ie))=1; in(mesh.edge(2,ie))=1; in=find(in);
    rect(ic,:)=((ismember(mesh.rect(1,:),in) & ismember(mesh.rect(2,:),in))  | ...
        (ismember(mesh.rect(2,:),in) & ismember(mesh.rect(3,:),in))  | ...
        (ismember(mesh.rect(3,:),in) & ismember(mesh.rect(4,:),in))  | ...
        (ismember(mesh.rect(4,:),in) & ismember(mesh.rect(1,:),in))) & ...
        ismember(mesh.rect(5,:),find(geom.semiconductor));
end
mesh.contact_nadj=(sum(rect,2)).';
mesh.contact_adj=zeros(geom.nc,max(sum(rect,2)));
for ic=1:geom.nc;
    if(mesh.contact_nadj(ic)>0), mesh.contact_adj(ic,1:sum(rect(ic,:)))=find(rect(ic,:)); end, end
%==============================================================================================100
% here we find the segments adjacent to contacts, used for schottky contacts
mesh.contact_seg=zeros(geom.nc,mesh.nn);
for ic=1:geom.nc;
    segment=find((ismember(geom.dgm(6,:),find(geom.reg2contact==ic))| ...
        ismember(geom.dgm(7,:),find(geom.reg2contact==ic)))& ...
        not(ismember(geom.dgm(6,:),0)|ismember(geom.dgm(7,:),0)));
    ie=ismember(mesh.edge(3,:),segment);
    i1=mesh.edge(1,ie); i2=mesh.edge(2,ie);
    x1=mesh.node(1,i1); x2=mesh.node(1,i2);
    y1=mesh.node(2,i1); y2=mesh.node(2,i2);
    temp=sqrt((x1-x2).^2+(y1-y2).^2);
    mesh.contact_seg(ic,:)=full(sparse(1,[i1 i2],[temp temp]/2,1,mesh.nn)); end
%**********************************************************************************************100
%
%
%**********************************************************************************************100
function [mesh]=f_dop(geom,mesh);
mesh.dop=zeros(1,mesh.nn);
mesh.dop_d=zeros(1,mesh.nn);
mesh.dop_a=zeros(1,mesh.nn);
%
for sd=1:geom.nd;
    if(geom.semiconductor(sd)),
        dgvet=geom.dgvet{sd};
        dtype=geom.dtype{sd};
        in1=mesh.rect(1,:); in2=mesh.rect(2,:);
        in3=mesh.rect(3,:); in4=mesh.rect(4,:);
        ir=(mesh.rect(5,:)==sd); ii=false(1,mesh.nn); ii([in1(ir) in2(ir) in3(ir) in4(ir)])=1;
        xrec=mesh.node(1,ii); yrec=mesh.node(2,ii);
        dop_node=f_EvalMolar(xrec,yrec,dgvet);
        %
        if strcmp(dtype,'N')
            mesh.dop_d(ii)=dop_node;
            mesh.dop_a(ii)=0;
        elseif strcmp(dtype,'P')
            mesh.dop_d(ii)=0;
            mesh.dop_a(ii)=dop_node;
        end
        mesh.dop(ii)=mesh.dop_d(ii)-mesh.dop_a(ii);
    end, end
%***********************************************************************************************100
%
%
%***********************************************************************************************100
function [mesh]=f_trap(geom,mesh);
ntrap=size(geom.dop_t,1);
mesh.dop_t=zeros(ntrap,mesh.nn);
for sd=1:geom.nd;
    if(geom.semiconductor(sd)),
        ir=find(mesh.rect(5,:)==sd);
        in1=mesh.rect(1,ir); in2=mesh.rect(2,ir);
        in3=mesh.rect(3,ir); in4=mesh.rect(4,ir);
        for ii=1:ntrap;
            dop_t=geom.dop_t(ii,sd);
            mesh.dop_t(ii,in1)=dop_t; mesh.dop_t(ii,in2)=dop_t;
            mesh.dop_t(ii,in3)=dop_t; mesh.dop_t(ii,in4)=dop_t;
        end, end, end
%***********************************************************************************************100
%
%
%***********************************************************************************************100
function [mesh]=f_ethj(geom,mesh);
% find edges on heterojunction segments
ie=ismember(mesh.edge(3,:),geom.ethj);
i1=mesh.edge(1,ie); i2=mesh.edge(2,ie);
% find nodes on heterojunction segments
inethj=false(1,mesh.nn);
inethj(i1)=1; inethj(i2)=1;
inethj=find(inethj);
nbrethj=length(inethj);
% add multiple points at heterointerfaces
mesh.node=[mesh.node mesh.node(:,inethj)];
% find rectangles with at least one node on the heterojunction segments
ir1=find(ismember(mesh.rect(1,:),inethj));
ir2=find(ismember(mesh.rect(2,:),inethj));
ir3=find(ismember(mesh.rect(3,:),inethj));
ir4=find(ismember(mesh.rect(4,:),inethj));
% now we must correct pointer to multiple nodes in the rect matrix
sd = geom.ethj_sd(2);
% regular node -> doubled node
temp=zeros(1,mesh.nn); temp(inethj)=(1:nbrethj)+mesh.nn;
if(mesh.rect(5,ir1)==sd), mesh.rect(1,ir1)=temp(mesh.rect(1,ir1)); end
if(mesh.rect(5,ir2)==sd), mesh.rect(2,ir2)=temp(mesh.rect(2,ir2)); end
if(mesh.rect(5,ir3)==sd), mesh.rect(3,ir3)=temp(mesh.rect(3,ir3)); end
if(mesh.rect(5,ir4)==sd), mesh.rect(4,ir4)=temp(mesh.rect(4,ir4)); end
% store pointers to multiple nodes
mesh.ethj=zeros(3,nbrethj);
mesh.ethj(1,:)=inethj;
mesh.ethj(2,:)=(1:nbrethj)+mesh.nn;
% update mesh.nn
mesh.nn=mesh.nn+nbrethj;
% store length of each segment along the heterointerface
x1=mesh.node(1,i1); x2=mesh.node(1,i2);
y1=mesh.node(2,i1); y2=mesh.node(2,i2);
temp=sqrt((x1-x2).^2+(y1-y2).^2);
ll=full(sparse(1,[i1 i2],[temp temp]/2,1,mesh.nn));
mesh.ethj(3,:)=ll(inethj);
%update mesh.edge!
%***********************************************************************************************100
%
%
%***********************************************************************************************100
function [mesh]=f_edge(geom,mesh);
nr=mesh.nr; % number of rectangles
edge=zeros(1,mesh.ne);
region=mesh.rect(5,:);
ir1=mesh.rect(10,:); ir2=mesh.rect(11,:); ir3=mesh.rect(12,:); ir4=mesh.rect(13,:);
ie1=mesh.rect(6,:); ie2=mesh.rect(7,:); ie3=mesh.rect(8,:); ie4=mesh.rect(9,:);
region1=zeros(1,nr); I1=not(ir1==0); region1(I1)=region(ir1(I1));
region2=zeros(1,nr); I2=not(ir2==0); region2(I2)=region(ir2(I2));
region3=zeros(1,nr); I3=not(ir3==0); region3(I3)=region(ir3(I3));
region4=zeros(1,nr); I4=not(ir4==0); region4(I4)=region(ir4(I4));
ir=1:nr;
in1=mesh.rect(1,:); in2=mesh.rect(2,:); in3=mesh.rect(3,:); in4=mesh.rect(4,:);
x=mesh.node(1,:); y=mesh.node(2,:);
x1=x(in1); x2=x(in2); x3=x(in3); x4=x(in4);
y1=y(in1); y2=y(in2); y3=y(in3); y4=y(in4);
xc=(x1+x2)./2; yc=(y2+y3)./2;
dl=geom.dl;
in1=edge; in2=edge;

for segment=1:size(geom.dgm,2);
    xs1=geom.dgm(2,segment); xs2=geom.dgm(3,segment);
    ys1=geom.dgm(4,segment); ys2=geom.dgm(5,segment);
    
    if(geom.dgm(7,segment)==0),
        sd1=geom.dgm(6,segment);
        sd2=geom.dgm(7,segment);
    else, sd1=geom.dgm(7,segment);
        sd2=geom.dgm(6,segment); end
    
    I=(region==sd1);
    
    iir=ir(I); r1=region1(I); r2=region2(I); r3=region3(I); r4=region4(I); r=region(I);
    xxc=xc(iir); yyc=yc(iir);
    
    if(not(dl(2)==dl(3)) & not(dl(4)==dl(5))), % obliquo
        I=(xxc>min(xs1,xs2) & xxc<max(xs1,xs2) & yyc>min(ys1,ys2) & yyc<max(ys1,ys2)); end
    
    if(dl(2,segment)==dl(3,segment)), % verticale
        I=(yyc>min([ys1 ys2]) & yyc<max([ys1 ys2]));
        iir=iir(I); r1=r1(I); r2=r2(I); r3=r3(I); r4=r4(I); xxc=xxc(I);
        [dmin,index]=min(abs(xxc-xs1));
        xmin=max([abs(x1(iir(index))-xs1) abs(x2(iir(index))-xs1)]);
        I=(abs(xxc-xs1)<xmin);
        iir=iir(I); r1=r1(I); r2=r2(I); r3=r3(I); r4=r4(I);
        I1=(r1==sd2); edge(ie1(iir(I1)))=segment;
        in1(ie1(iir(I1)))=mesh.rect(3,iir(I1));
        in2(ie1(iir(I1)))=mesh.rect(2,iir(I1));
        I3=(r3==sd2); edge(ie3(iir(I3)))=segment;
        in1(ie3(iir(I3)))=mesh.rect(1,iir(I3));
        in2(ie3(iir(I3)))=mesh.rect(4,iir(I3)); end
    
    if(dl(4,segment)==dl(5,segment)), % orizzontale
        I=(xxc>min([xs1 xs2]) & xxc<max([xs1 xs2]));
        iir=iir(I); r1=r1(I); r2=r2(I); r3=r3(I); r4=r4(I); yyc=yyc(I);
        [dmin,index]=min(abs(yyc-ys1));
        ymin=max([abs(y1(iir(index))-ys1) abs(y4(iir(index))-ys1)]);
        I=(abs(yyc-ys1)<ymin);
        iir=iir(I); r1=r1(I); r2=r2(I); r3=r3(I); r4=r4(I);
        I2=(r2==sd2); edge(ie2(iir(I2)))=segment;
        in1(ie2(iir(I2)))=mesh.rect(2,iir(I2));
        in2(ie2(iir(I2)))=mesh.rect(1,iir(I2));
        I4=(r4==sd2); edge(ie4(iir(I4)))=segment;
        in1(ie4(iir(I4)))=mesh.rect(4,iir(I4));
        in2(ie4(iir(I4)))=mesh.rect(3,iir(I4)); end, end

ne=mesh.ne; mesh.edge=zeros(4,ne);
mesh.edge(1,:)=in1; mesh.edge(2,:)=in2; mesh.edge(3,:)=edge;

% apply boundary conditions on edges
for ii=1:size(geom.dgm,2);
    I=find(mesh.edge(3,:)==ii);
    mesh.edge(4,I)=geom.edge(ii); end
%***********************************************************************************************100
%
%
%***********************************************************************************************100
function [mesh]=rectmesh2(geom);
% node matrix
div_x=ones(1,length(find(not(cellfun('isempty',geom.X))))-1);
div_y=ones(1,length(find(not(cellfun('isempty',geom.Y))))-1);
I=find(not(cellfun('isempty',geom.X))); J=find(not(cellfun('isempty',geom.Y)));
xb=evalgeom(geom,geom.X(I)); yb=evalgeom(geom,geom.Y(J));
Nx=sum(div_x)+1; Ny=sum(div_y)+1; nrx=Nx-1; nry=Ny-1; nr=nrx*nry; nn=Nx*Ny;
x=zeros(1,Nx); y=zeros(1,Ny);
index=1;
for k=1:length(div_x),
    x(index:sum(div_x(1:k))+1)=linspace(xb(k),xb(k+1),div_x(k)+1);
    index=index+div_x(k);
end
index=1;
for k=1:length(div_y),
    y(index:sum(div_y(1:k))+1)=linspace(yb(k),yb(k+1),div_y(k)+1);
    index=index+div_y(k);
end
mesh.xgrid=x; mesh.ygrid=y;
[x,y]=meshgrid(x,y);
ix=1:Nx; iy=1:Ny; [ix,iy]=meshgrid(ix,iy);
in=(ix-1)*Ny+iy; % nodal matrix index
in1=reshape(in(2:Ny,1:Nx-1)',1,nr);
in2=reshape(in(2:Ny,2:Nx)',1,nr);
in3=reshape(in(1:Ny-1,2:Nx)',1,nr);
in4=reshape(in(1:Ny-1,1:Nx-1)',1,nr);
x1=x(in1); x2=x(in2); x3=x(in3); x4=x(in4);
y1=y(in1); y2=y(in2); y3=y(in3); y4=y(in4);
l1=y2-y3; l2=x2-x1; l3=y1-y4; l4=x3-x4;
xc=(x1+x2)./2; yc=(y2+y3)./2;
x=reshape(x,1,nn); y=reshape(y,1,nn);
mesh.node=[x; y];
% rect matrix
region=zeros(1,nr);
nreg=max([geom.dgm(6,:) geom.dgm(7,:)]);
for sdl=1:nreg;
    [xv,yv]=extractpatch(geom.dgm,sdl);
    I=inpolygon(xc,yc,xv,yv); I=find(I); region(I)=sdl;
end
iemat=reshape(cumsum(ones(1,(Ny-1)*Nx)),Nx,Ny-1)';
ie1=reshape(iemat(:,2:Nx)',1,nr);
ie3=reshape(iemat(:,1:Nx-1)',1,nr);
iemat=reshape(cumsum(ones(1,Ny*(Nx-1))),Nx-1,Ny)';
iemat=iemat+(Ny-1)*Nx;
ie4=reshape(iemat(1:Ny-1,:)',1,nr);
ie2=reshape(iemat(2:Ny,:)',1,nr);
irmat=zeros(Ny+1,Nx+1);
irmat(2:Ny,2:Nx)=reshape(cumsum(ones(1,nr)),Nx-1,Ny-1)';
ir1=reshape(irmat(2:Ny,3:Nx+1)',1,nr);
ir3=reshape(irmat(2:Ny,1:Nx-1)',1,nr);
ir4=reshape(irmat(1:Ny-1,2:Nx)',1,nr);
ir2=reshape(irmat(3:Ny+1,2:Nx)',1,nr);
mesh.rect=[in1; in2; in3; in4; region; ie1; ie2; ie3; ie4; ir1; ir2; ir3; ir4];
mesh.nr=nr; mesh.nrx=nrx; mesh.nry=nry;
mesh.nn=nn; mesh.nnx=Nx; mesh.nny=Ny;
mesh.ne=(Nx-1)*Ny+Nx*(Ny-1);
% edge matrix
mesh=f_edge(geom,mesh);
%***********************************************************************************************100
