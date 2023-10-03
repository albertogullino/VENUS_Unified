tp=1;
for k=1:730
 tp=Tstor(:,:,k)*tp;
 if k>702
  map(tp)
  'k = ' , k
  pausak
 end
end 

