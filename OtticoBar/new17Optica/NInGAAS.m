function nindex=nInGaAs(lambda,x)
%NALGAAS(LAMBDA,X) Returns the refractive index of AlGaAs
%       LAMBDA is the wavelength in m (can be a vector) and
%       X is the compositional Al content.
%       NALGAAS(LAMBDA,X) returns an array with the same size as LAMBDA
%
%       Original model from M.A. Afromovitz, Solid State Comm, 15, 1974
%       Model extended with the range 0.8µm - egamma from
%       S.W. Corzine, Ph.D. Dissertation

%if min(lambda) < 800e-9
    l=lambda;
    h=6.6261*10^-34;
    c=2.9979*10^8;
    q=1.6022*10^-19;

    e0=1.43-1.53*x+0.45*x.^2;
    de0=0.341-0.09*x+0.45*x.^2;
    
    
    ed0=e0+de0;

    cc=h*c/(l)/q./e0;
    cs0=h*c/(l)/q./ed0;

    a0=6.3-1.16*x;
    b0=9.4+0.75*x;

    m=cc;
    fcc=m.^(-2).*(2-sqrt(1+m)-sqrt(1-m));
    m=cs0;
    fcs0=m.^(-2).*(2-sqrt(1+m)-sqrt(1-m));
    eps=real(a0.*(fcc+1/2*(e0./ed0).^(3/2).*fcs0)+b0);
    nindex=(sqrt(eps));


