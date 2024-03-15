function [y] = HF_diagH(qD,theta)
%
global PPEc PPEv
% global PPVcc PPVvv % non-screened potentials
global PPWcc PPWvv % screened potentials
global K1 EFc EFv kBTev C0 V0 flagLimit
%
K2 = K1 + qD.*exp(1i*theta);
k2 = abs(K2);
LqD=log(qD);
%
%Ek2c = ppval(PPEc,k2);
Ek2v = ppval(PPEv,k2);
%
%fi=find(Ek2c>=abs(C0));
fpes1=ones(size(k2));
fpes2=fpes1;
%if flagLimit==1
%    fpes1(fi)=0;
%end
fi=find(Ek2v>=abs(V0));
if flagLimit==1
    fpes2(fi)=0;
end

%rhok2c = fpes1./(1.0D0 + exp((Ek2c-EFc)/kBTev));
rhok2v = fpes2./(1.0D0 + exp((Ek2v-EFv)/kBTev));
%
% y =  (ppval(PPVcc,qD).*rhok2c + ...
%       ppval(PPVvv,qD).*rhok2v ) .* qD; % non-screened potentials
%Wc=-exp(ppval(PPWcc,LqD));
Wv=-exp(ppval(PPWvv,LqD));
%
%y =  (Wc.*rhok2c + Wv.*rhok2v )  .* qD; % screened potentials
y =  ( Wv.*rhok2v )  .* qD; % screened potentials
