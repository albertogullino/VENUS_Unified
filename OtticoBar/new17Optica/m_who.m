DR=whos;
l=length(DR);
for kD=1:l
 dun=getfield(DR(kD),'size');
 if max(dun)>100
  pausak
 end

end
