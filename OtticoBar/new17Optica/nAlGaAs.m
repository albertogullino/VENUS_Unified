function nindex=nAlGaAs(lambda,x)
%NALGAAS(LAMBDA,X) Returns the refractive index of AlGaAs
%       LAMBDA is the wavelength in m (can be a vector) and
%       X is the compositional Al content.
%       NALGAAS(LAMBDA,X) returns an array with the same size as LAMBDA
%
%       Original model from M.A. Afromovitz, Solid State Comm, 15, 1974
%       Model extended with the range 0.8µm - egamma from
%       S.W. Corzine, Ph.D. Dissertation

%'cont', keyboard
if size(x,1)>1
 x=x';
end
%'cont', keyboard
%if min(lambda) < 800e-9
if min(lambda) < 740e-9
   disp('Too short wavelength.');
   return
elseif max(lambda) > 1.2e-6
   disp('Too large wavelength.');
   return
end

 if lambda < 800e-9
% if lambda < 900e-9
    l=lambda;
    h=6.6261*10^-34;
    c=2.9979*10^8;
    q=1.6022*10^-19;

    e0=1.425+1.155*x+0.37*x.^2;
    ed0=1.765+1.115*x+0.37*x.^2;

    cc=h*c/(l)/q./e0;
    cs0=h*c/(l)/q./ed0;

    a0=6.3+19.0*x;
    b0=9.4-10.2*x;

    m=cc;
    fcc=m.^(-2).*(2-sqrt(1+m)-sqrt(1-m));
    m=cs0;
    fcs0=m.^(-2).*(2-sqrt(1+m)-sqrt(1-m));
    eps=a0.*(fcc+1/2*(e0./ed0).^(3/2).*fcs0)+b0;
    nindex=sqrt(eps);

 else

nindex=zeros(max(size(x)),max(size(lambda)));
e1=1.239852195651227e-06./lambda;
e0=3.65+.871*x+.179*x.^2;
ed=36.1-2.45*x;
egamma=1.424+1.266*x+.26*x.^2;
egammam = egamma - 25e-3;
ef=(2*e0.^2-egamma.^2).^(.5);
eta=pi*ed.*(2*(e0.^3).*(e0.^2-egamma.^2)).^(-1);
m1=eta.*(2*pi)^(-1).*(ef.^4-egamma.^4);
m3=eta/pi.*(ef.^2-egamma.^2);

indexLower=find(e1>egammam);
indexHigher=find(e1<=egammam);

%'qui', keyboard

nindex(indexHigher)=(1+m1(indexHigher)+...
   m3(indexHigher).*e1.^2+...
   eta(indexHigher)/pi.*(e1.^4).*...
   log((ef(indexHigher).^2-e1.^2)./...
   (egamma(indexHigher).^2-e1.^2))).^(.5);

ngammam(indexLower)=(1+m1(indexLower)+...
   m3(indexLower).*egammam(indexLower).^2+...
   eta(indexLower)/pi.*(egammam(indexLower).^4)...
   .*log((ef(indexLower).^2-egammam(indexLower).^2)./...
   (egamma(indexLower).^2-egammam(indexLower).^2))).^(.5);
nindex(indexLower)=ngammam(indexLower)+.35*(e1*ones(size(indexLower))-...
   egammam(indexLower));

 end

nindex=real(nindex);
if x==0
%'index', keyboard
end
