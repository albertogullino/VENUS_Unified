%
%function is=ismat(M)
%
% is=0  M scalar
% is=1  M vector
% is=-1 M row
% is=n  M matrix n dim

function is=ismat(M)

s=size(M);
if min(s)==1
 is=1;
 if s(1)==1
  is=-1;
 end
 if max(s)==1
  is=0;
 end
else
 is=length(s);
end

