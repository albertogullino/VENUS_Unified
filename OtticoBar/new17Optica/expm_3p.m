imetex=-1;
imetex=2;
if imetex==0
x=expm(x);
elseif imetex==1
[x,D] = eig(x);
x = x * diag(exp(diag(D))) / x;
%'ver d', keyboard
elseif imetex==2
 P=x;
 TscattuuM
elseif imetex==-1
 Oo=x; 
end