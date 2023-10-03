 function [f1,zi,nz,az]=fieldz_siyi(dx,r,KKi,f1,zi,nz,Nx,az,a)
 if nargout==1
  ize=0;
 else
  ize=1;
end
 ic=length(f1(1,:));
 for nxi=1:Nx
  ic=ic+1;
  if ize==1
   nz(ic,:)=r;
   az(ic,:)=a;
   zi=[zi dx];
  end 
  f1(:,ic)=KKi*f1(:,ic-1);
 end
