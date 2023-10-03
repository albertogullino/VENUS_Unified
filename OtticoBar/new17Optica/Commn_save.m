function []=Common_save

'entro comm;', keyboard
GG=whos('global');

stringa = 'save(''glo.mat'','
stringaG = 'global '

for ind = 1:length(GG)
    stringa = [stringa, '''',GG(ind).name,''','];
    stringaG = [stringaG, GG(ind).name,' '];
end

stringa=[stringa, '''stringaG''',');'];

'verifico'
keyboard


eval(stringaG)
eval(stringa)
