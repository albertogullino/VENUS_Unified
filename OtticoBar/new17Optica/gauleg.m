function  [xi,wi]=gauleg(X1,X2,N)

    EPS=1.D-7;
    PIG=pi;

    M=(N+1)/2;
    XM=(X1+X2)/2;
    XL=(X2-X1)/2;

    for I=1:M
     Z=cos(PIG*(I-.25)/(N+.5));
     Z1=2*Z;
     while abs(Z-Z1)>EPS
        P1=1;
        P2=0;
       for J=1:N
        P3=P2;
        P2=P1;
        P1=((2.D0*J-1.D0)*Z*P2-(J-1.D0)*P3)/J;
       end  %J
       PP=N*(Z*P1-P2)/(Z*Z-1);
       Z1=Z;
       Z=Z1-P1/PP;
     end  %while
     xi(I)=XM-XL*Z;
     xi(N+1-I)=XM+XL*Z;
     wi(I)=2*XL/((1-Z*Z)*PP*PP);
     wi(N+1-I)=wi(I);
    end   %I
