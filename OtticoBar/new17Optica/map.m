%function h=map(M,x,y,ifi,ibar)

function h=map(M,x,y,ifi,ibar)

%'entro map', keyboard

if nargin>=4
 figure(ifi)
else
 h=figure;
end
  COL='colorbar';
if nargin==5
 if ibar==1
  COL='colorbar(''horiz'')';
 elseif ibar==0
  COL='colorbar(''vert'')';  
 else 
  COL='colorbar off';
 end
end
M=real(M);
warning off
if nargin>1
  s=size(M);
  sx=size(x);
  sy=size(y);
 if length(sx)~=length(sy)
   ' error xy sizes ' ,
   keyboard
   return
 elseif min(sx)==1
  if length(x)==s(2)
   surf(x,y,M), shading('interp'), view(0,90), eval(COL)
  else
   surf(y,x,M), shading('interp'), view(0,90), eval(COL)
  end
 else
  if sx==s
   surf(x,y,M), shading('interp'), view(0,90), eval(COL)
  else
   ' error xy sizes ' ,
   keyboard
   return
  end
 end
else
   surf(M), shading('interp'), view(0,90), eval(COL)
end
axis square

%'fine map', keyboard