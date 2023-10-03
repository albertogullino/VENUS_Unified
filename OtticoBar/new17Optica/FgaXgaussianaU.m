function [Gdbr]=gaDBR(kx,ky,pg);

lambda=pg.lambda;
ring=pg.r_in;
rout=pg.r_out;
Nb=pg.ND;


nvh=2.3;
Lvh=lambda/4/nvh;

nvl=1.3;
Lvl=lambda/4/nvl;

Lb=[Lvl Lvh];
nb=[nvl nvh];

la=lambda;
rr=1;

kr=2*pi/lambda*rr;
kv=sqrt(kx^2+ky^2);
kt=kv/kr;
%pausak
if kt>.9
 kt=.9*kr;
 %' kt > 1', 
 %' kt > 1', keyboard
end

[Ge,Gm]=gaperdm(kt,0,la,Lvh,nvh,Lb,nb,Nb,rout,rr,0,[],[],ring);

Gdbr(1,1)=Ge;
Gdbr(2,2)=Gm;