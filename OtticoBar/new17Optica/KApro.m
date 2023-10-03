 [vv,ee]=eig(P);
 ei=j*imag(diag(ee));
 Od=vv*diag(exp(ei))*inv(vv);
 Mn=Od;
 
 I=eye(size(Ga1));
 D1=(prod(diag(I-Ga1.^2)))^(.25/length(KK));
 D2=(prod(diag(I-Ga2.^2)))^(.25/length(KK));
 V1=[I Ga1; Ga1  I]/D1; 
 V2=[I -Ga2; -Ga2  I]/D2; 
 MnR=V2*Mn*V1;
