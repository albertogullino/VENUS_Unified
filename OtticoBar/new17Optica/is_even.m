%
% function is=is_even(n)
%
% is=1 if n is even
%
function is=is_even(n)
for k=1:length(n)
isf=n(k)/2-fix(n(k)/2);
if isf==0
 is(k)=1;
else
 is(k)=0;
end
end
