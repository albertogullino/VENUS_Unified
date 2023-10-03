%clear 
%close all
%load grat
' entro le_nogra'

 fi=find(th>0);
 tha=th;
 nla=nl;
 
   [dilua,ailua,nilua,fstua]=lens_add(tha,h,r,N,nla,ifps,Npair);
if ifp==-10
   figure, plot(ailua,cumsum(dilua),'r.'), pausak
end   
   
    fi=find(th<0);
    thb=abs(th(fi));
    nlb=[nl(1) nl(fi+1) nl(end)];
    
 if Np_ad~=0 
  [dilub,ailub,nilub,fstub]=lens_add(thb,h,r,N,nlb,ifps,Np_ad);
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

%' fine no return ', keyboard
%figure, plot(cumsum(dilua),real(nilua(:,1)),'r.',cumsum(dilub(fi)),real(nilub(fi,1)),'g.'), pausak
%figure, plot(ailub(fi,:),csa+cumsum(dilub(fi)),'w.',ailua,cumsum(dilua),'.'), pausak

sa=size(ailua,2);
sb=size(ailub,2);
if sa~=sb
 if sa>sb
  le=length(dilub);
  ze=zeros(le,1);
  nilub=[nilub ze];
  ailub=[ailub ze];
 else
  le=length(dilua);
  ze=zeros(le,1);
  nilua=[nilua ze];
  ailua=[ailua ze];
 end
end

dilu=[dilua; dilub(fi)];
ailu=[ailua; ailub(fi,:)];
nilu=[nilua; nilub(fi,:)];
fstu=[fstua; fstub(fi,:)];



if ifp==-10
   fim=find(fstu(:,2)>1);
   df=diff(fim);
   if length(find(df>1))>0
    nrep=fstu(fim([1 end]),2);
   else
    nrep=fstu(fim(1),2);
   end 


fi=find(dilu>1e-5);
 adisp=ailu(fi,:);
 ddisp=dilu(fi);
 ndisp=nilu(fi,:);
 fdisp=fstu(fi,:);
 
    for k=1:length(nrep)
     fim=find(fdisp(:,2)==nrep(k));
     fumu=fim;
     nd=ndisp(fim,:);
     ad=adisp(fim,:);
     dd=ddisp(fim);
     fd=fdisp(fim,:);
     nd0=nd;
     ad0=ad;
     dd0=dd;
     fd0=fd;
     for kk=1:nrep(k)-1
      fumu=[fumu; fim+length(fim)*kk];
      nd=[nd; nd0];
      ad=[ad; ad0];
      dd=[dd; dd0];
      fd=[fd; fd0];
     end
     fisp1=1:fim(1)-1;
     fisp2=fim(end)+1:length(ddisp);
     ndisp=[ndisp(fisp1,:); nd; ndisp(fisp2,:)];
     adisp=[adisp(fisp1,:); ad; adisp(fisp2,:)];
     ddisp=[ddisp(fisp1); dd; ddisp(fisp2)];
     fdisp=[fdisp(fisp1,:); fd; fdisp(fisp2,:)];

%     adisp
%     adisp(fim
    end 

   csu1=cumsum(ddisp);
   figure, plot(adisp,csu1,'r.'),
   hold on
   fim=find(fdisp(:,2)>1);
   nrep=fdisp(fim([1 end]),2);
   
   col='cg';
   for k=1:length(nrep)
     fim=find(fdisp(:,2)==nrep(k));
     plot(adisp(fim,:),csu1(fim),[col(k),'o']),
   end
   title([' repetition ',num2str(nrep')])   
   pausak


% per vedere il profilo di indice   
  n0=ndisp(:,1);
  fc=find(diff([1; n0])~=0);
  fc=[fc; length(n0)];
  dd0=csu1(1:fc(1)-1);
  nn0=n0(1:fc(1)-1);
  for k=1:length(fc)-1
     dd0=[dd0; csu1(fc(k)-1); csu1(fc(k):fc(k+1)-1)];
     nn0=[nn0; n0(fc(k)); n0(fc(k):fc(k+1)-1)];
  end
  %   dd0=[dd0; csu1(fc(k):length(n0))];
  %   nn0=[nn0; n0(fc(k):length(n0))];

  
    figure,  plot(dd0,nn0,csu1(fc-1),n0(fc-1),'.')
    pausak
    figure, plot(cumsum(dilu),nilu,'.')
     pausak
end