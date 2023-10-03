clear
close all
load pri

lf=length(phiv);
al=phiv'*pi/180;
df=diff(al(1:2))/pi;

mv=[0  2  4];
mv=[0  2  4]+1;
mv=[1 3 5]-1;
%mv=[1 ];



ii=4;

for ii=1:4
ii
pausak

im=0;
clear T11
for m=mv
 im=im+1;
 imp=0;
 for mp=mv
  imp=imp+1;
  ff=reshape(Gtem1(:,1,ii),lf,1);
  fa=repmat(-sin(al*m).*sin(al*mp),1,1);
  fa([1 end])=  fa([1 end])/2;
  fint=sum(ff.*fa)*df;
  Tem(im,imp)=fint;
  
  ff=reshape(Gtm1(:,1,ii),lf,1);
  fa=repmat(cos(al*m).*sin(al*mp),1,1);
  fa([1 end])=  fa([1 end])/2;
  fint=sum(ff.*fa)*df;
  %T11(2,2)=fint;  
  Tmm(im,imp)=fint;
  
  ff=reshape(Gte1(:,1,ii),lf,1);
  fa=repmat(-sin(al*m).*cos(al*mp),1,1);
  fa([1 end])=  fa([1 end])/2;
  fint=sum(ff.*fa)*df;
  Tee(im,imp)=fint;
  
  ff=reshape(Gtme1(:,1,ii),lf,1);
  fa=repmat(cos(al*mp).*cos(al*m),1,1);
  fa([1 end])=  fa([1 end])/2;
  fint=sum(ff.*fa)*df;
  Tme(im,imp)=fint;

%  figure, plot(phiv,ff.*fa), pausak
 end
 
end 

Tt11=[Tee Tem; Tme Tmm];

sum(sum(abs(Tt11)))
pausak
end


