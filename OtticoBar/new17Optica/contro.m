clear
close all
load pri

lf=length(phiv);
al=phiv'*pi/180;
df=diff(al(1:2))/pi;

mv=[0  2  4];
mv=[0  2  4]+1;
mv=1;
im=0;

for m=mv
 im=im+1;
 imp=0;
 for mp=mv
  imp=imp+1;
  ff=reshape(Gte1(:,1,1:4),lf,4);;
  fa=repmat(sin(al*m).*sin(al*mp),1,4);
  fa=repmat(cos(al*m).*cos(al*mp),1,4);
%  fa=repmat(sin(al*m).*cos(al*mp),1,4);
%  fa=repmat(cos(al*m).*sin(al*mp),1,4);
  fint=sum(ff.*fa)*df;
  figure, plot(phiv,ff.*fa), pausak
  %fint=sum(reshape(Gtm1(:,1,1:4),lf,4).*repmat(sin(al*m).*cos(al*mp),1,4))*df;
  %fint=sum(reshape(Gtem1(:,1,1:4),lf,4).*repmat(sin(al*m).*cos(al*mp),1,4))*df;
  %fint=sum(reshape(Gtme1(:,1,1:4),lf,4).*repmat(sin(al*m).*cos(al*mp),1,4))*df;
  PI1(im,imp)=fint(1);
  PI2(im,imp)=fint(2);
  PI3(im,imp)=fint(3);
  PI4(im,imp)=fint(4);
 end
end 
fi=find(abs(PI1)<1e-10); 
PI1(fi)=0;
fi=find(abs(PI2)<1e-10); 
PI2(fi)=0;
fi=find(abs(PI3)<1e-10); 
PI3(fi)=0;
fi=find(abs(PI4)<1e-10); 
PI4(fi)=0;

real(PI1), 'PI1', pausak
real(PI2), 'PI2', pausak
real(PI3), 'PI3', pausak
real(PI4), 'PI4', pausak

PT=[PI1  PI3; PI4 PI2]