idue=0;
if exist('gra_le')
 if isfield(gra_le,'pos')==1
 if gra_le.pos<=2
  idue=1;
 end
 end
end 


if idue==1
 n0s=nv0;
 
 nsos=gra_le.n_pa;
 nv0=n0s;
 fi_zgrat
 la_pa=lambda;

 nsos=gra_le.n_ve;
 nv0=n0s;
 fi_zgrat
 la_ve=lambda; 
 
 lambda=max([la_pa la_ve]);
 lambdavet=sort([la_pa la_ve]);

 
 
 
else
'qui fiz', keyboard
 fi_z
 lambda_ve=lambda;
 la_ve=lambda;
 la_pa=lambda;
 
end