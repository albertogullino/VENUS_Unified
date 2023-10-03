function []=Common_save

GG=whos('global');

stringa = 'save(''glo.mat'',''GG'',' ;
stringaG = 'global ';

for ind = 1:length(GG)
    stringa = [stringa, '''',GG(ind).name,''','];
    stringaG = [stringaG, GG(ind).name,' '];
end

stringa=[stringa, '''stringaG''',');'];
%'ve', keyboard
eval(stringaG)
eval(stringa)
