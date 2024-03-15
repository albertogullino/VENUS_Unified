function [epsfun] = f_EvalLindhardScreening(kgrid,wk,wth,VK,VTH,vic,viv)
%
global PPEc PPVcc PPEv PPVvv
global PPEcV PPVccV PPEvV PPVvvV
global qD NQW
%
epsfun = zeros(1,length(kgrid));
%
hLindhard = @(k,theta) Lindhard_c(k,theta);
%
for ic = vic;
    PPVcc = PPVccV{ic};
    PPEc = PPEcV{ic};
    %
    for iq = 1:length(kgrid)
        qD = kgrid(iq);
        hLindhardIntegrand=hLindhard(VK,VTH);
        hLindhardIntegral=2*wk*hLindhardIntegrand*wth';
        epsfun(iq) = epsfun(iq) + 1/(2*pi)^2 * hLindhardIntegral;
    end
end
%
hLindhard = @(k,theta) Lindhard_v(k,theta);
for iv = viv
    PPVvv = PPVvvV{iv};
    PPEv = PPEvV{iv};
    for iq = 1:length(kgrid)
        qD = kgrid(iq);
        hLindhardIntegrand=hLindhard(VK,VTH);
        hLindhardIntegral=2*wk*hLindhardIntegrand*wth';
        epsfun(iq) = epsfun(iq) + 1/(2*pi)^2 * hLindhardIntegral;
    end
end
%
epsfun = 1 - epsfun/NQW;
%
%==========================================================================
function [y] = Lindhard_v(k,theta)
%
global PPEv PPVvv
global EFv
global V0
global flagLimitMaxDensity
global kBTev qD qel
%
K = k.*exp(1i*theta);
%
Ek1v = ppval(PPEv,k);
Ek2v = ppval(PPEv,abs(K-qD));

fi=find(Ek1v>=abs(V0));
fpes1=ones(size(k));
fpes2=fpes1;
if flagLimitMaxDensity==1
    fpes1(fi)=0;
end
fi=find(Ek2v>=abs(V0));
if flagLimitMaxDensity==1
    fpes2(fi)=0;
end
%
rhok1v = fpes1./(1 + exp((Ek1v-EFv)/(kBTev)));
rhok2v = fpes2./(1 + exp((Ek2v-EFv)/(kBTev)));
%
y = ppval(PPVvv,qD) .* (rhok2v - rhok1v)./(qel*(Ek2v - Ek1v)) .* k;
%==========================================================================
function [y] = Lindhard_c(k,theta)
%
global PPEc PPVcc
global EFc
global C0
global kBTev qD qel
global flagLimitMaxDensity
%
K = k.*exp(1i*theta);
%
Ek1c = ppval(PPEc,k);
Ek2c = ppval(PPEc,abs(K-qD));

fi=find(Ek1c>=abs(C0));
fpes1=ones(size(k));
fpes2=fpes1;
if flagLimitMaxDensity==1
    fpes1(fi)=0;
end
fi=find(Ek2c>=abs(C0));
if flagLimitMaxDensity==1
    fpes2(fi)=0;
end
%
rhok1c = fpes1./(1 + exp((Ek1c-EFc)/(kBTev)));
rhok2c = fpes2./(1 + exp((Ek2c-EFc)/(kBTev)));
%
y = ppval(PPVcc,qD) .* (rhok2c - rhok1c)./(qel*(Ek2c - Ek1c)) .* k;
%==========================================================================