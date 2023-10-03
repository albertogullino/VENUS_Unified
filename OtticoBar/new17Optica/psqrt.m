function w=psqrt(z)
%determination appropriate for passive media: imag{kz}<0
w=sqrt(z);
i=find(imag(w)>0);
w(i)=-w(i);
return