function [vet,val]=eig_mio(A,B,nsol);
%' Eig mio in'
%pausak

%   [vet,vali]=eig(B,-A);
%   val=diag(vali);
%   return

global pMu0u pMu0u1 lKA nk1max nures pMc sim0 numodi pMei iLP ldap

if nargin==2
 nsol=0;
end

iveme=0;

if nsol>0
 iveme=2;
end

%A=N1;
%B=N2;
vet=zeros(length(A),nsol);

%Kplot1=A;
%Kplot2=B;
%an_mat1
iaccp=1;
npias=100;

nn0=length(find(iaccp>0));
siac=size(iaccp);
if siac(1)==0
 nn0=npias;
 siac=1;
 iaccp=1:nn0;
 iaccv=1:nn0;
% iveme=2;
else
 nn0=length(find(iaccp>0));
% iveme=0;
end

if iveme==0
 siac=1;
end

%' iacc eig', keyboard

val=[];
vet=[];
ics=0;

for ks=1:siac(1)

 if length(iaccp(ks,:))==npias | iveme==0
  N1=A;
  N2=B;
  puv=1:length(A);
  nso=length(A);
 else
    puu=[];
    puv=[];
    ad0=1:nk1max;
    fi=find(iaccp(ks,:)~=0);
    iaccv=iaccp(ks,fi);
    if length(iaccv)>0
     for ipu=iaccv
      adv=(ipu-1)*nk1max;
      puv=[puv ad0+adv];
      for ipui=iaccv
       add=((ipu-1)*lKA+ipui-1)*nk1max;
       puu=[puu pMei+add];
      end
     end
     puu=sort(puu);
     sim=sim0*length(iaccv)/2;
     N1=reshape(A(puu),sim,sim);
     N2=reshape(B(puu),sim,sim);
    end
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

 val=[val; diag(vali)];
 vet(puv,ics+[1:nso])=veti;
 ics=ics+nso;

end  %ks

%disp('eig mio'), keyboard
fi0=find(abs(val)>0);
val=val(fi0);
vet=vet(:,fi0);

[du,is]=sort(1./abs(val));
val=val(is);
vet=vet(:,is);
%' Eig mio out'
%pausak
