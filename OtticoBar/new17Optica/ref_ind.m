%
%
% function n=ref_ind(ifmat,lambda,x,y)
%
% Ref. inded of materials
%
% ifmat=1: Al_x Ga_1-x As
% ifmat=2: In_x Ga_1-x As
%
% ifmat and x with the same size
%
% lambda in meter

function n=ref_ind(nvv,ifmat,lambda,x,y)

n=zeros(size(ifmat))';
nvv0=nvv;
for ima=-4:2
 fim1=find(ifmat==ima);
 if length(fim1)>0
  if ima==1
   n(fim1)=nAlGaAs(lambda,x(fim1));
  elseif ima==2
   n(fim1)=ningaas(lambda,x(fim1));
  elseif ima==-1
   [nIRe,keIm]=nAu(lambda);
   nvv0=0;
   n(fim1)=nIRe-i*keIm;
  elseif ima==-2
   [nIRe,keIm]=nTi(lambda);
   n(fim1)=nIRe-i*keIm;
   nvv0=0;
  elseif ima==-3
   [nIRe,keIm]=nPt(lambda);
   n(fim1)=nIRe-i*keIm;
   nvv0=0;
  elseif ima==-4
   [nIRe,keIm]=nCr(lambda);
   n(fim1)=nIRe-i*keIm;
   nvv0=0;
  end
 end
end
n=n.'+i*imag(nvv0);
