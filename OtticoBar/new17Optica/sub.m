global ireset_int
icont_it=1;
if length(ireset_int)==0
 ireset_int=-1;
end
%'quo', keyboard
if ireset_int==-1
% sub_old
 sub_oldpeak
else
 sub_int
end
