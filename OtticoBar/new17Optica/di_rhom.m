
nv=4;
fia=2*pi/nv;
nfia=501;
xt=linspace(0,2*pi,nfia);
dae=-0.0;
dap=-0.06;
 ru=1+dap*cos(nv*xt)+dae*cos(3*nv*xt);
 figure, polar(xt,ru,'r'),

