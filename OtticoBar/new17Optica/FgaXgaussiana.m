function [Gdbr]=gaDBR(kx,ky,lambda);



nvh=2.3;
Lvh=lambda/4/nvh;

nvl=1.3;
Lvl=lambda/4/nvl;

Lb=[Lvl Lvh];
nb=[nvl nvh];

la=lambda;
Nb=11;
ring=1;
rr=1;
rout=1;

kr=2*pi/lambda*rr;
kv=sqrt(kx^2+ky^2);
kt=kv/kr;
%pausak
if kt>1
 kt=kr;
 %' kt > 1', 
 %' kt > 1', keyboard
end

[Ge,Gm]=gaperdm(kt,0,la,Lvh,nvh,Lb,nb,Nb,rout,rr,0,[],[],ring);

Gdbr(1,1)=Ge;
Gdbr(2,2)=Gm;