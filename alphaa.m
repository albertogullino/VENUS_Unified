alpha=1e-2; % um-1

Lz=linspace(0,100,100);
f=exp(-alpha.*Lz);

figure,plot(Lz,f)

R=mode.Rrad;
MM=[Se1.*R(in1) Se2.*R(in2) Se3.*R(in3)];
PspBulk=1000*sum(sparse(MM))*h*(Clight*1e-2/mean(1e-9*mode.vlambda))

Rrad2D=mode.RradQW;
MM = [Lp1.*Rrad2D(iiQW1) Lp2.*Rrad2D(iiQW2)];
frsp=mode.frsp;
Psp=frsp*1000*sum(sparse(MM))*h*(Clight*1e-2/mean(1e-9*mode.vlambda))