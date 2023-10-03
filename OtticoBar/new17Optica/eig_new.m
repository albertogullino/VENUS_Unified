function [vet,val]=eig_mio(A,B,nsol);
%' Eig mio in'
%pausak


global pMu0u pMu0u1 lKA nk1max nures pMc sim0 numodi pMei iLP
global PUp1 PUp2 polsub

if nargin==2
 nsol=0;
end

iveme=0;

if nsol>0
 iveme=2;
end

%' new', keyboard
if polsub==-1
 if length(PUp1)==0
  Kplot1=A;
  an_mat1
  PUp1=PUm;
 else
  PUm=PUp1;
 end
else
 if length(PUp2)==0
  Kplot1=A;
  an_mat1
  PUp2=PUm;
 else
  PUm=PUp2;
 end
end



val=zeros(size(A));
vet=val;

for ks=1:length(PUm)

 if length(PUm)==1
  N1=A;
  N2=B;
  pum=1:length(A)^2;
 else
   pum=PUm{ks};
   sim=sqrt(length(pum));
   N1=reshape(A(pum),sim,sim);
   N2=reshape(B(pum),sim,sim);
 end

  if iveme>=1
   N1i=-inv(N1);
   Mtot=N1i*N2;

   if iveme==1
    [veti,vali]=eig(Mtot);
   elseif iveme==2
    OPTIONS.disp=0;
    nso=fix(length(iaccv)/nn0*nsol)
    if nso<2
     nso=4;
    end
    [veti,vali]=eigs(Mtot,nso,OPTIONS);
%    [veti,vali]=eigs(Mtot,nso,'LR',OPTIONS);
%    [vett,valt]=eig(Mtot);
   end

  elseif iveme==0   %iveme

   [veti,vali]=eig(N2,-N1);

  end
 val(fix(pum))=vali;
 vet(fix(pum))=veti;

end  %ks

%disp('eig mio'), keyboard
val=diag(val);
fi0=find(abs(val)>0);
val=val(fi0);
vet=vet(:,fi0);
[du,is]=sort(1./abs(val));
val=val(is);
vet=vet(:,is);
%' Eig mio out'
%pausak
