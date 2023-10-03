function W=F2(A,Mvefm0,Mvefp0,Mvegm0,Mvegp0,besp,besm,xdx1,xdx3)


   AnFs=A;

   Mvefm=diag(AnFs)* Mvefm0;
   Mvegm=diag(AnFs)*Mvegm0;

   Mvefp=diag(AnFs)*Mvefp0;
   Mvegp=diag(AnFs)*Mvegp0;
   
   Exdu=besp'*Mvefp+besm'*Mvefm;
   Eydu=besp'*Mvegp+besm'*Mvegm;
   
   I=mean(abs(Exdu.^2+Eydu.^2),2);
   Es=xdx3*I;
   En=xdx1*I;
   W=sqrt(2*Es./En);


