function ierro=ver_par(ch_p)
s=length(ch_p);
ierro=0;
for ke=1:s
 ste=ch_p{ke,1};
 for ki=ke+1:s
  if isequal(ste,ch_p{ki,1})==1
   ' errore parametri '
   ierro=1;
  end
 end
end
