function screen(string,number,string2,number2,string3,number3)

if nargin==2
 disp([string,'  ' num2str(number)])
elseif nargin==4
 disp([string,'  ' num2str(number), string2,'  ' num2str(number2)])
elseif nargin==6
 disp([string,'  ' num2str(number), string2,'  ' num2str(number2), string3,'  ' num2str(number3)])
end

