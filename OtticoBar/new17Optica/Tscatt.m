%      [xve,Dau] = eig(xx);
[xve,Dau] = eig(xx);
ve=xve;

diaut=diag(Dau);
if abs(real(diaut(1)))>abs(imag(diaut(1)))
 [du,fis]=sort(real(diaut));
 diaut=diaut(fis);
 ve=ve(:,fis);
end

vei=inv(ve);
lve=length(ve)/2;
l1=1:lve;
l2=l1+lve;

analitico =diaut(l1);

Ta=vei;
Tb=Ta;
Ta11=Ta(l1,l1);
Ta12=Ta(l1,l2);
Ta21=Ta(l2,l1);
Ta22=Ta(l2,l2);
iT=inv(Ta22);
sb11=-iT*Ta21;
sb12=iT;
sb21=Ta11-Ta12*iT*Ta21;
sb22=Ta12*iT;

Ta=ve;

Ta11=Ta(l1,l1);
Ta12=Ta(l1,l2);
Ta21=Ta(l2,l1);
Ta22=Ta(l2,l2);

iT=inv(Ta22);
sa11=-iT*Ta21;
sa12=iT;
sa21=Ta11-Ta12*iT*Ta21;
sa22=Ta12*iT;


uno=diag(ones(lve,1));
De=inv(uno-sa22*sb11*exp(2*analitico));
s21u=sa21*sb21*exp(analitico)*De;
s12u=sa12*sb12*exp(analitico)*De;
s11u=sa11+sb11*sa12*sa21*exp(2*analitico)*De;
s22u=sb22+sa22*sb12*sb21*exp(2*analitico)*De;

is=inv(s12u);
Tu(l1,l1)=s21u-s22u*is*s11u;
Tu(l1,l2)=s22u*is;
Tu(l2,l1)=-is*s11u;
Tu(l2,l2)=is;
xax=Tu.';