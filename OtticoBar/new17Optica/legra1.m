clear 
close all
load grat

Np0=Npair;
%tha0=[tha(1:2)];
%nla0=[nla(1:3) nla(2)];
[dilua0,ailua0,nilua0,fstua0]=lens_add(thgsu*5,tha,h,r,N,nla,ifps,Np0);

   fi=find(fstua0(:,2)==1);   
   csu=cumsum(dilua0);
   figure, plot(ailua0,csu,'.',ailua0(fi),csu(fi),'r.'), pausak

tha1=[tha];
nla1=[nla(3) nla(2:end)];
Np1=Npair+2;
   [dilua1,ailua1,nilua1,fstua1]=lens_add(thgsu,tha1,h,r,N,nla1,ifps,Np1);
   
   fi=find(fstua1(:,2)==1);   
   csu=cumsum(dilua1);
   figure, plot(ailua1,csu,'.',ailua1(fi),csu(fi),'r.'), pausak
   
   return
   di0=dilua(fi);
   ai0=ailua(fi);
      cs0=cumsum(di0);

   figure, plot(ai0,cs0,'c.'), 
