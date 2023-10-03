function [E,F]=f_EvalRefractiveIndexExpansions(NHarmonics,csi,d,alpha1,a1,x1,n1,n2)
%
% [E,F]=f_EvalRefractiveIndexExpansions(NHarmonics,csi,d,alpha1,a1,x1,n1,n2)
%
%------------------------------------------------------------------------
% Fourier/Floquet transform of n^2(x), 1/n^2(x)
%------------------------------------------------------------------------
% In this function the coefficients of the Fourier/Floquet transforms of
% n^2(x) and 1/n^2(x) are computed, being n^2(x) a raised-cosine function.
% The degeneration in rect case is accounted for.
%
% NHarmonics is number of Floquet harmonics, csi the related kx, d is the
% period of the structure, alpha1 the roll-off factor of the raised-cosine
% (alpha1=0 means rect), a1 the width, x1 the displacement from x=0, n1,n2
% the levels.
%
% Alberto Tibaldi, 04/03/2015
%------------------------------------------------------------------------

for indr=1:NHarmonics
    for indc=1:NHarmonics
        vs=csi(indr)-csi(indc);
        if(vs==0)
            E(indr,indc)=n2^2+2/d*a1*(n1^2-n2^2); 
            F(indr,indc)=1./(n2^2)+2/d*a1*(1./(n1^2)-1./(n2^2));
        elseif(abs(abs(2*alpha1*vs*a1/pi)-1)<=1e-10)
            E(indr,indc)=2/d*a1*sin(vs*a1)/(vs*a1)*(n1^2-n2^2)*(pi/4)*exp(j*vs*x1);
            F(indr,indc)=2/d*a1*sin(vs*a1)/(vs*a1)*(1./(n1^2)-1./(n2^2))*(pi/4)*exp(j*vs*x1);
         else
            E(indr,indc)=2/d*a1*sin(vs*a1)/(vs*a1)*(n1^2-n2^2)*cos(alpha1*vs*a1)/(1-(2*alpha1/pi*vs*a1)^2)*exp(j*vs*x1);
            F(indr,indc)=2/d*a1*sin(vs*a1)/(vs*a1)*(1./(n1^2)-1./(n2^2))*cos(alpha1*vs*a1)/(1-(2*alpha1/pi*vs*a1)^2)*exp(j*vs*x1);
        end
    end
end

return