% QW plots
x=mesh.xgrid*1e4;
xQW=x(1:mesh.nnxQW{1});
y=mesh.ygrid*1e4;

figure,plot(xQW,mesh.T(mesh.inMQW{1}),xQW,mesh.T(mesh.inMQW{3}))
title('Temp in QWs'),xlabel('\rho, \mum'),ylabel('\DeltaT, K')
pausak
figure,plot(xQW,mode.nQW{end}{1},xQW,mode.nQW{end}{2},xQW,mode.nQW{end}{3})
title('Electrons in QWs'),xlabel('\rho, \mum'),ylabel('Carriers level, cm^{-2}')
legend('QW1','QW2','QW3')
pausak
figure,plot(xQW,mode.pQW{end}{1},xQW,mode.pQW{end}{2},xQW,mode.pQW{end}{3})
title('Holes in QWs'),xlabel('\rho, \mum'),ylabel('Carriers level, cm^{-2}')
pausak
Rc=reshape(mode.ecb,mesh.nny,mesh.nnx);
Rv=reshape(mode.evb,mesh.nny,mesh.nnx);

figure, plot(y,Rv(:,[1 mesh.nnxQW{1}])), hold on, plot(y,Rc(:,[1 mesh.nnxQW{1}]))
title('Energy band diagram, \rho=0 - \rho=\rho_{QW}'),xlabel('z, \mum'),ylabel('Energy, eV')
xlim([110 119])

n2D=zeros(mesh.nnxQW{1},length(mode.ii_dd));
p2D=n2D;

for i2D=1:length(mode.ii_dd)
    n2D(:,i2D)=MODEplot{end}.nQW{i2D}{2};
    p2D(:,i2D)=MODEplot{end}.pQW{i2D}{2};
end

WQW=mesh.vWMQW{1};

figure
subplot(121)
hold on, grid on,box on
plot(xQW,n2D/WQW,'.-','linewidth',2)
plot(xQW,p2D/WQW,'.-','linewidth',2)
xlabel('\rho, \mum'),ylabel('QW carrier densities, cm^{-2}')
subplot(122)
hold on, grid on,box on
plot(mode.ii_dd*1e3,n2D(1,:)/WQW,'.-','linewidth',2)
plot(mode.ii_dd*1e3,p2D(1,:)/WQW,'.-','linewidth',2)
xlabel('Current, mA'),ylabel('QW carrier densities, cm^{-3}')
legend('n','p','location','northwest')

El2=reshape(mode.N2D,mesh.nny,mesh.nnx);
Ho2=reshape(mode.P2D,mesh.nny,mesh.nnx);
figure, semilogy(y,El2(:,1),y,Ho2(:,1))


El=reshape(mode.elec,mesh.nny,mesh.nnx);
Ho=reshape(mode.hole,mesh.nny,mesh.nnx);
figure, semilogy(y,El2(:,1)+El(:,1),y,Ho2(:,1)+Ho(:,1))
title('Bulk+Bound carrier densities at \rho=0')
xlabel('z, \mum'),ylabel('Carrier Densities, cm^{-3}')
axis([114 116 1e10 1e20]),grid on

rox=2.175
% ir=sum(geom.div_x(1:3));
[~,ir]=min(abs(x-rox));

[~,fiy]=min(abs(y-mesh.yMQW{2}*1e4));
[~,fix]=min(abs(x-xQW(end)));

figure, plot(xQW,El(fiy,1:fix),xQW,Ho(fiy,1:fix)),
grid on,hold on
title('Bulk carrier densities at central QW')
xlabel('\rho, \mum'),ylabel('Carrier densities, cm^{-3}')
pausak
 figure,plot(ir*mesh.node(2,mesh.ICAV),'ro')
 zmin=min(mesh.node(2,mesh.ICAV));
 zmax=max(mesh.node(2,mesh.ICAV));
 fi=find(mesh.ygrid>=zmin & mesh.ygrid<=zmax );
 figure, plot(y), hold on, plot(fi,y(fi),'r.') 


% tag punti cavità
nnxCAV = mesh.nnxQW{1};
nnyCAV = length(mesh.ICAV)/nnxCAV;
ICAV = reshape(mesh.ICAV,nnyCAV,nnxCAV);
figure,plot(mesh.node(1,mesh.inMQW{1}),mode.elec(ICAV)) % lungo rho
figure,plot(mesh.node(1,mesh.inMQW{1}),mode.hole(ICAV)) % lungo rho
figure,plot(mesh.node(2,mesh.inMQW{1}),mode.elec(ICAV),'.') % lungo z


ficav=fi;
figure, plot(x,El(ficav,:)), pausak
figure, plot(x,Ho(ficav,:)), pausak


JN_X=reshape(mode.Jn_x,mesh.nny,mesh.nnx);
JN_Y=reshape(mode.Jn_y,mesh.nny,mesh.nnx);
JP_X=reshape(mode.Jp_x,mesh.nny,mesh.nnx);
JP_Y=reshape(mode.Jp_y,mesh.nny,mesh.nnx);

NQW=mesh.NMQW;

RRaug=(reshape(mode.RAuger,mesh.nny,mesh.nnx));
RRsrh=(reshape(mode.RSRH,mesh.nny,mesh.nnx));
RRrad=(reshape(mode.Rrad,mesh.nny,mesh.nnx));
Cnd=(reshape(mode.Ccapn3D,mesh.nny,mesh.nnx))/NQW*1e-9*mesh.vWMQW{1};
Cpd=(reshape(mode.Ccapp3D,mesh.nny,mesh.nnx))/NQW*1e-9*mesh.vWMQW{1};

Rtot=RRaug+RRsrh+RRrad;
Rtn=Rtot-Cnd;
Rtp=Rtot-Cpd;

J_X=JN_X+JP_X;
J_Y=JN_Y+JP_Y;
J_Xn=JN_X;
J_Yn=JN_Y;

leQW=length(mode.nQW{end}{2});
xQW=x(1:leQW);
    
figure, plot(x,-J_Xn(fiy,:),'linewidth',2),
hold on
plot(xQW,-mode.JnQW{2},'linewidth',2), 
xlabel('radial coord.')
ylabel('Electron flux')
legend(' 3D ',' 2D')
pausak
%figure, plot(y,Rc-Rv)
dX=100e-7;
xlim([mesh.yMQW{3}-dX,mesh.yMQW{1}+dX]*1e4)
title(' Bands in Cavity')

figure,hold on,plot(x,VELMInfo(end).E2), pausak

% plots doping e Free Carrier
% figure,plot(mesh.ygrid,mode.dop_dp(1:mesh.nny),mesh.ygrid,mesh.dop_d(1:mesh.nny),mesh.ygrid,mode.elec(1:mesh.nny)+10,'k'),set(gca,'yscale','log')
 %figure,plot(mesh.ygrid,mode.dop_dp(1:mesh.nny),mesh.ygrid,mesh.dop_d(1:mesh.nny),mesh.ygrid,mode.elec(1:mesh.nny)+10,'k'),legend('ND+','ND','elec'),set(gca,'yscale','log')
% figure,plot(mesh.ygrid,mode.dop_dp(1:mesh.nny),mesh.ygrid,mesh.dop_d(1:mesh.nny),mesh.ygrid,mode.elec(1:mesh.nny)+10,'k'),set(gca,'yscale','log'),ylim([1e17,1e20])
 
 figure,plot(mesh.ygrid,mode.dop_dp(1:mesh.nny),mesh.ygrid,mesh.dop_d(1:mesh.nny),mesh.ygrid,mode.elec(1:mesh.nny)+10,'k'),
 legend('ND+','ND','elec','location','best'),set(gca,'yscale','log')
 figure,plot(mesh.ygrid,mode.dop_am(1:mesh.nny),mesh.ygrid,mesh.dop_a(1:mesh.nny),mesh.ygrid,mode.hole(1:mesh.nny)+10,'k'),
 legend('AM+','AM','hole','location','best'),set(gca,'yscale','log')
 
  figure,ind=25;plot(mesh.ygrid,mode.ecb(mesh.nny*ind+[1:mesh.nny]),'b',mesh.ygrid,mode.evb(mesh.nny*ind+[1:mesh.nny]),'r',mesh.ygrid,mode.EFn(mesh.nny*ind+[1:mesh.nny]),'k--',mesh.ygrid,mode.EFp(mesh.nny*ind+[1:mesh.nny]),'k-.','LineWidth',2)


% potenziale

PHI=reshape(mode.phi,mesh.nny,mesh.nnx);
 figure, contourf(mesh.X*1e4,mesh.Y*1e4,PHI)