%  r_out
%       ____    n4   ____  ^
%       | n2|  n1 |     |      | thick
% __ _|     |___ |     |      v
%         d2   d1
%            n3
% r_in

function [Te,Tm,Glate,Glatm,Eze,Ezm,Hze,Hzm,nef]=Orta_tracar(tetai,r_in,r_out,r1,r2,d1i,d2i,thicki,lambdai,Nmodi,itetm,Strut,ifp);


%if nargin<13
 ifield=0;
%end

ico=0;

if nargin<11
 isto=0;
end
if nargin<12
 itetm=3;
end
d1=d1i*1000;
d2=d2i*1000;

lambda=lambdai*1000;
thick=thicki*1000;

theta=tetai*pi/180;

n1=r1;
n2=r2;
%n1=r_out;
n3 = real(r_in);
n4 = real(r_out);

n3 = r_in;
n4 = r_out;

d=d1+d2;


Nstratper=length(thick); % numero di strati in ogni dente


if length(thick)~=length(n2)
 'errore parametri in orta_gen: thick non consistente con n1,n2',  keyboard
end


errore = 1e-6;
formulaz='HFIE';

mu0 = 4e-7*pi; 
c = 2.997925e8;   %nm / ns, quindi frequenze in GHz

if itetm==1 | itetm==3

%' prima di suborta', keyboard
polariz='TE';
sub_orta
%neff(1)=bet(1)/k0;
%Glate= s11;
if ifield==1
 Eze=Ez;
 Hze=Hz;
else
 Eze=0;
 Hze=0;
end

Te=Tr;
Glate=s11;
nef(1)=bet(1)*1e-3/k0;
%'beta TE', keyboard

else
Glate=0;
Te=0;
Eze=0;
Hze=0;

end

if itetm==2 | itetm==3

polariz='TM';
sub_orta

Tm=Tr;
Glatm=s11;
if ifield==1
 Ezm=Ez;
 Hzm=Hz;
else
 Ezm=0;
 Hzm=0;
end
nef(2)=bet(1)*1e-3/k0;

%'beta TM', keyboard

else
Tm=0;
Glatm= 0;
Ezm=0;
Hzm=0;
end


if ifp==-10
' dentro orta_tracar campi', keyboard
end
