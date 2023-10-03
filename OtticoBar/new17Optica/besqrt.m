function f=besqrt(x)

f=sqrt(x);
fip=find(imag(f)>0);
if length(fip)>0
 f(fip)=conj(f(fip));
end
