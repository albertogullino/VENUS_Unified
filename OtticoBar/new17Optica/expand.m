%
%function  res=expand(ns,rep,nsnot)
%
% ns :   vettore di numero strati da ripetere
% rep:   vettore di ripetizione
% nsnot: valori da escludere nella ripetizione
%
% res: puntatore espanso

function  res=expand(ns,rep,nsnot)

 l=length(rep);
 cr=1;
 k=1;
 while k<=l
  if rep(k)==1
   res(cr,1)=k;
   cr=cr+1;
   k=k+1;
  else
   if length(find(nsnot==rep(k)))==0
    rei=rep(k);
   else
    rei=1;
   end
   nsi=ns(k);
   for krep=1:rei
    res(cr:cr+nsi,1)=[k:k+nsi]';
    cr=cr+nsi;
   end
   k=k+nsi;
  end
 end
if res(end)>length(ns)
 res=res(1:end-1);
end
