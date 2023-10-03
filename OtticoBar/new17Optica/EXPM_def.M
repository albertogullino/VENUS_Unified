function E=expm_def(M)
[x,D] = eig(M);
E = x * diag(exp(diag(D))) / x;
