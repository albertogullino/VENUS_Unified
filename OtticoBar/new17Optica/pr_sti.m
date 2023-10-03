clear
close all
d=.2;
D=1;
ns=21;
X=10;
Y=10;
x=linspace(-X,X,100);
y=linspace(-1,Y,100);
figure
xv=[-X X X -X];
yv0=[-d -d d d]/2;
for k=1:ns
yv=yv0+(k-1-fix(ns/2))*D;
%xv
%yv
%pausak
 fill(xv,yv,'w'), hold on
end

npr=1001;
rv=linspace(0,2*D+d/2,npr);
rv=linspace(0,X,npr);

hold on, polar(linspace(0,2*pi,100),max(rv)*ones(1,100))
axis equal
pausak


nubesu=10;

muv=[0:2:2*nubesu];
AB=ones(length(rv),length(muv))*NaN;
im=0;
for mu=muv

im=im+1;

 if mu==0

  fat=pi/2;
  F(1)=fat;
  for k=2:npr
   r=rv(k);
   a0=d/(2*r);
   nstr=fix((r+d/2)/D);
   som=0;
   if nstr==0
    if a0<=1
     as0=asin(a0);
    else
     as0=fat;
    end
    som0=as0;
   else
    som0=asin(a0);
   end
   for kj=1:nstr
    a1=(kj*D+d/2)/r;
    if a1<=1
     as1=asin(a1);
    else
     as1=fat;
    end
    a2=(kj*D-d/2)/r;
    as2=asin(a2);
    som=som+(as1-as2);
   end
   F(k)=som+som0;
  end  %r
  ver=4*(F*rv')*diff(rv(1:2))/(pi*max(rv)^2)/(d/D)
  pausak

 else
  fat=0;
  F(1)=fat;
  for k=2:npr
   r=rv(k);
   a0=d/(2*r);
   nstr=fix((r+d/2)/D);
   som=0;
   if nstr==0
    if a0<=1
     as0=asin(a0);
    else
     as0=fat;
    end
    som0=sin(mu*as0)/mu;
   else
    as0=asin(a0);
    som0=sin(mu*as0)/mu;
   end

   for kj=1:nstr
    a1=(kj*D+d/2)/r;
    if a1<=1
     as1=asin(a1);
    else
     as1=fat;
    end
    a2=(kj*D-d/2)/r;
    as2=asin(a2);
    som=som+(sin(mu*as1)-sin(mu*as2))/mu;
   end
   F(k)=som+som0;
  end  %r
  F=F;
 end
 AB(:,im)=4*F';

% figure, plot(rv,F,'.g')
% pausak
end
 figure, plot(rv,AB)
