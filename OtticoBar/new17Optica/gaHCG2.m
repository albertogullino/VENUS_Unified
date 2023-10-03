function [Ga,Tr]=gaHCG2(kx,ky,pg);

lambda=pg.lambda;
DC=pg.DC;
d=pg.period;
d1i=d*DC;
d2i=d-d1i;
thick=pg.thick;
%NModi=11;
NModi=pg.NModi;
r1=pg.n1;
r2=pg.n2;
r_in=pg.r_in;
r_out=pg.r_out;


kr=sqrt(kx^2+ky^2);
k0=2*pi/lambda;
if kx>.9*k0
 %'ke>k0', keyboard
 kx=.9*k0;
end
if ky>.9*k0
 %'ke>k0', keyboard
 ky=.9*k0;
end
%tetai=asin(kr/k0);
%phii=atan2(ky,kx);

%[tetai phii], pausak

 %'ke>k0', keyboard
  [Ga,Tr]=orta_skewTr(kx,ky,r_in,r_out,r1,r2,d1i,d2i,thick,lambda,NModi,0);
  %[Ge,Gm]=bast_genS(tetai,r_in,r_out,r1,r2,d1i,d2i,thick,lambda,itetm,Nmodi);  
%  Gtm20(ife,ite)=Gtmd(2);
%  Geo(ite)=Glate;
%  Gmo(ite)=Glatm;
%  pausak

