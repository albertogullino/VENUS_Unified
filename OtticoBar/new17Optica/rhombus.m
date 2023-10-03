fi=linspace(0,2*pi,100);
mo=ones(size(fi));
dap=.05;
a=5.2/2;
b=5.7/2;

          nv_sh=4;
          xt=fi;
          rapax=b/a;
          dae=(rapax-1)/2;
          ru=(1+dap*cos(nv_sh*xt)+dae*cos(2*xt));
          roc2=mo.*ru.*exp(j*fi);

rou=(roc2-(min(real(roc2))+j*min(imag(roc2))))*a;
dx=max(real(rou));
dy=max(imag(rou));
disp(' size(x,y) ')
[dx dy]
figure, plot(rou), axis equal, grid
