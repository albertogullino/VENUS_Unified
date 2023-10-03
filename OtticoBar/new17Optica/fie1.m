load sa1
Nx=50;
%Nx=2;

[fiez,nz,zet,Gaqw,NQW_ef,az]=Flam_field(L_i,n_i1,iat,fiQW,la,rr,kt,ibast,par_grat,ga,Nx,a_i);

Ez=sum(fiez);
Ez=Ez/max(Ez);
gatot=ga;
ga=gatot/NQW_ef;

zqw=sum(L_i(1:fiQW(1)))/1000;
[du,imi]=min(abs(zet-zqw));
I0=abs(Ez).^2;
I=log10(I0/I0(imi))+3;
figure, plot(zet,I, zet, nz(:,1),'r'), 
 title(['  lambda_{res} = ',num2str(la),' Gth= ',num2str(ga)]), 
  xlabel(' Long. coord. (um)')
 ylabel(' Intens (Log10) / Index')
 a=axis;
 a(3)=-3;
 a(4)=5;
 axis(a)
 pausak
 figure, plot(zet,I0*3.52, zet, nz(:,1),'r'), 
  title(['  lambda_{res} = ',num2str(la),' Gth= ',num2str(ga)]), 
   xlabel(' Long. coord. (um)')
  ylabel(' Intens  / Index')
