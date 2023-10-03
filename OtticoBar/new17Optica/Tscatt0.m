%      [xve,Dau] = eig(xx);
[xve,Dau] = eig(xx);
ve=xve;
diaut=diag(Dau);
fze=1;
if abs(real(diaut(1)))>1
% 'passo 0 aut'
 fze=0;
end
vei=inv(ve);
xax=ve*diag(exp(diaut))*vei*fze;
