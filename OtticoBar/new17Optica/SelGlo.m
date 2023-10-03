%function []=Common_save


stringaGi=stringaG;
nul='&';

for ind = 1:length(GG0)
 gi=[GG0(ind).name,' '];
 fis=findstr(stringaGi,gi);
 fiN=fis+[0:length(gi)-1];
 stringaGi(fiN)=repmat(nul,1,length(fiN));
 
 %length(stringaGi)
 %pausak
end

for k=1:length(stringaGi)
 if strcmp(stringaGi(k),'&')
  stringaGi(k)=' ';
 end
end

stringaclear=['clear ',stringaGi];
save saclearglo stringaclear

%load saclearglo
%eval(stringaclear)
%'ve', keyboard
