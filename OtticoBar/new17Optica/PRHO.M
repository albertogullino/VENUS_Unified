bvero=3;
avero=2;
dap=.04;
dap=.1;

nv=4;
fia=2*pi/nv;
nfia=201;
xt=linspace(0,2*pi,nfia);

rapax=bvero/avero;
dae=-(1-rapax)/(1+rapax)*(1+dap);
R00=bvero/(1+dae+dap);
ru=R00*(1+dap*cos(nv*xt)+dae*cos(2*xt));

Rvd=ru.*exp(j*xt);
figure, plot(Rvd), grid, axis equal
