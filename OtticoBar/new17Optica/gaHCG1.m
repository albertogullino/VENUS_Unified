function [Ga,Tr]=gaHCG1(kx,ky);

lambda=1.616;
DC=.7;
d=.7;
d1i=d*DC;
d2i=d-d1i;
thick=.35;
NModi=11;
r1=3.48;
r2=1;
r_in=1;
r_out=1;
itetm=1;  %1 TE, 2 TM, 3 entrambe

kr=sqrt(kx^2+ky^2);
k0=2*pi/lambda;
if kr>k0
 %'ke>k0', keyboard
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

