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
xax=ve*diag(exp(diaut))*vei;
if max(max(abs(xax)))>100
 xax=zeros(size(xax));
end
