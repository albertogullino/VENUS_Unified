load stopLoop
%istopLoop
%istopLoop=input(' input istop Loop = ? ');
if istopLoop==1
 istopLoop=0
else
 istopLoop=1
end
if length(istopLoop)==1
  save stopLoop istopLoop -v4
end

