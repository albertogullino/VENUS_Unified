%
%
% function [mue,muh,me,mh,tc]=momatc(ifmat,T,x,y)
%
% Mobility of electrons and holes (mue, muh)
% Masses of electrons and holes (me, mh)
% Thermal conductivity
%
% ifmat=1: Al_x Ga_1-x As
% ifmat=2: In_x Ga_1-x As
%
% ifmat and x with the same size
%
%

function [mue,muh,me,mh,tc]=momatc(ifmat,T,x,y)

for ima=1:1
 fim1=find(ifmat==ima);
 if length(fim1)>0
  if ima==1
   [mue(fim1,1),muh(fim1,1),me(fim1,1),mh(fim1,1),tc(fim1,1)]=mmtalga(x(fim1),T);
%  elseif ima==2
%   n(fim1)=ningaas(lambda,x);
  end
 end
end

