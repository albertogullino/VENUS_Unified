function PTherm = HeatIntegral(mesh,HeatSource)

rho=mesh.xgrid*1e4;
drho=diff([0 rho]);
xdxN=rho.*drho;

z=mesh.ygrid*1e4;
dzN=diff([0 z(2:end-1)]);

HeatSource=reshape(HeatSource,mesh.nny,mesh.nnx);
qtot=HeatSource(2:end-1,:);
qtot(isnan(qtot))=0;

PTherm=2*pi*xdxN*qtot'*dzN'*1000; % mW