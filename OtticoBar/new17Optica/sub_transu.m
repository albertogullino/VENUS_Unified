%'su', keyboard
Ldum=L_i;
[du,fimal]=max(L_i);
Ldum(fimal)=max(dto);
 cs=cumsum(Ldum);
if ~exist('iwo')
 iwo=0;
end
if exist('ires_scatt')
 if ires_scatt==0
   if iwo==0
    subca_hcg  
   else 
%    subca_hcgw  
    subca_SUE 
   end
  else
   sub_transunuu   %metodo scattering
 end
else
 subca_hcg  
end


%subca_hcgu  