%
%krv=[.2 .4]
%tolku=.05;
Pus=[];
Pusf=[];
ldap=[0];

% if length(find(nmv==nm))==0 | ik~=1
%' proar ', keyboard
for nm=1:npam:numodi
 for ik=1:length(krv)
  kin=krv(ik)-tolku;
  kfi=krv(ik)+tolku;
%  kin=k_ret;
%  kfi=alim;
%  'modifica'
%  keyboard
%  if nm==1 & ik==1
  if length(find(nmv==nm))==1 & ik==1
%   k_red=[0 alim];
   k_red=0;
  else
   k_red=k_ret;
  end
%  k_red
%  keyboard
  pui=pu+(nm-1)*lk;
  Du=abs(Dm(pui));
  Du=Du/max(Du);
  der=diff([0; Du]);
 if length(k_red)==0 | k_red==0
   if ikr==3
    if nm==1
     fiu=pu;
    else
     fiu=find(KK>=kin);
    end
   else
%    fiu=find(KK<=kfi);
    fiu=find((KK<=alim) | (KK>=kin & KK<=kfi));
   end
   pua=fiu;
    fiuf=find((KK<=alim));
   puaf=fiuf;
%   'rimo', keyboard
 else
   puaf=[];
%  if nm==1 & ik==1
  if length(find(nmv==nm))>0 & ik==1
   fiu=find(KK>k_red & Du-pre>0);
  else
   fiu=find(KK>kin & Du-pre>0 & KK<=kfi);
  end
   if length(fiu)>4
    pua=fiu;
   else
    pua=lk-3:lk;
   end
  end
% if nm~=1 | ik~=1 | (ik==1 & nm==1)
  if npam==1
    Pusf=[Pusf pui(puaf)];
    Pus=[Pus pui(pua)];
    ldap=[ldap ldap(end)+length(pua)];
  else
   Pusf=[Pusf pui(puaf) pui(puaf)+lk];
   Pus=[Pus pui(pua) pui(pua)+lk];
   ldad=ldap(end)+length(pua);
   ldap=[ldap ldad ldad+length(pua)];
  end
%  ' pus', keyboard
% end

 end  %krv

end
