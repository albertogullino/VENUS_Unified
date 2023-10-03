clear 
close all
load grat

 'primo ', keyboard
   [dilua,ailua,nilua,fstua]=lens_add(thgsu,tha,h,r,N,nla,ifps,Npair);
    fi=find(th<0);
    thb=abs(th(fi));
    nlb=[nl(1) nl(fi+1) nl(end)];
  'secondo ', keyboard
 if Np_ad~=0 
  [dilub,ailub,nilub,fstub]=lens_add(thgsb,thb,h,r,N,nlb,ifps,Np_ad);
 else
  dilub=[];
  ailub=[];
  nilub=[];
  fstub=[];
 end
 %  [dilu,ailu,nilu,fstu]=lens_sav(th,h,r,N,nl,ifp,Hre,Npair,Rel,Rm_ring);
 

%figure, plot(ailub(fia,:),cumsum(dilub(fia)),'.',ailub(fis,:),0.2+cumsum(dilub(fis)),'.',ailub(fia,:),cumsum(dilub(fia)),'.',ailub(fiu,:),hr+cumsum(dilub(fiu)),'.'),
csa=sum(dilua);
%fi=find(cumsum(dilub)>h);
%fito=fi;
%fi=[fi(1)-1; fi];
if Np_ad~=0
 fidu=find(diff(fstub(:,2))~=0);
 fi=fidu(1)+1:length(fstub);
else 
 fi=[];
end 

' fine no return ', keyboard
figure, plot(cumsum(dilua),real(nilua(:,1)),'r.',csa+cumsum(dilub(fi)),real(nilub(fi,1)),'g.'), 
hold on
plot(cumsum(dilua),real(nilua(:,2)),'m.',csa+cumsum(dilub(fi)),real(nilub(fi,2)),'c.'), 

pausak
figure, plot(ailub(fi,:),csa+cumsum(dilub(fi)),'w.',ailua,cumsum(dilua),'.'), pausak