%      [xve,Dau] = eig(xx);
[xve,Dau] = eig(xx);
ve=xve;

diaut=diag(Dau);
if abs(real(diaut(1)))>abs(imag(diaut(1)))
 [du,fis]=sort(real(diaut));
 diaut=diaut(fis);
 ve=ve(:,fis);
 ' ass ',
else
 [du,fis]=sort(imag(diaut));
 diaut=diaut(fis);
 ve=ve(:,fis);
 ' prop ',
end

vei=inv(ve);
lve=length(ve)/2;
l1=1:lve;
l2=l1+lve;

analitico =diaut(l1);

Ta=vei;
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

E1=diag(exp(analitico));
E2=diag(exp(2*analitico));

uno=diag(ones(lve,1));
De1=inv(uno-sb11*sa22*E2);
De2=inv(uno-sa22*E2*sb11);

s11u=sa11+sa12*De1*sb11*sa21;
s12u=sa12*De1*sb12*E1;
s21u=sb21*E1*De2*sa21;
s22u=sb22+sb21*E1*De2*sa22;

is=inv(s12u);
Tu(l1,l1)=s21u-s22u*is*s11u;
Tu(l1,l2)=s22u*is;
Tu(l2,l1)=-is*s11u;
Tu(l2,l2)=is;
xax=Tu.';