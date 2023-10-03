%[x,D] = eig(x);
%x = x * diag(exp(diag(D))) / x;
x=expm(x);
%'ver d', keyboard
