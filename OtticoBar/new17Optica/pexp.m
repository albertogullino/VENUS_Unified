A=[1 j*.8; -1 pi];
[V,D0] = eig(A);
E1 = V * diag(exp(diag(D0))) / V;
f=sqrt(pi);
A=[1 j*.8; -1 pi]*f;
[V,D] = eig(A);
E2 = V * diag(exp(diag(D))) / V;
Et=E2*E1;
Etv = V * diag(exp(diag(D0*(1+f)))) / V;



