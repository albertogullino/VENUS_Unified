 load comparison

%   gpla=2e4*pi*rr/lambda*imag(Ksi)/Gam; 

   x=(hz')*1e-3;
   rperm=real(sqrt(relPerm));
    figure, plot(x,rperm,'r',x,Fi*4,'w')
    title([' lambda_{res} = ',num2str(lambda),'  Gth = ',num2str(gpla/3)]), pausak