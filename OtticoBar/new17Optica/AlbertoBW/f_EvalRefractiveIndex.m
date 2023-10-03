function [n]=f_EvalRefractiveIndex(x,d1,n1,n2,alpha1,x1)

n=sqrt(n2^2+(n1^2-n2^2)*f_cosrialz(x,x1,d1,alpha1));

return