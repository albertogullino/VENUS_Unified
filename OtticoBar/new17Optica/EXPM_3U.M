[Oo,D] = eig(Oo);
Oo = Oo * diag(exp(diag(D))) / Oo;
