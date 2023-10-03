load sa1
 yal=csus+sh;
   fia=find(ailuc(:,1)>0);
   pul=fia(end);
   yi=csu(pul)-h;
   yu=yal(end);
   ytot=csu;
   yu=ytot(end);
   [du,fi0]=min(abs(ytot-yi));
   yyy=ytot(fi0);
fi0s=fi0;
   clear atot
   atot=ailuc(1:fi0,:);
   ntot=niluc(1:fi0,:);
   ftot=fstuc(1:fi0,:);
%fi0=fi0-1;
   while  yyy<yu
    fi0=fi0+1;
    yyy=ytot(fi0);
%    [du,fial]=min(abs(yyy-yal));
    dif=yal-yyy;
    fis=find(dif>0);
    fial=fis(1);
  %  'in', pausak

%    ytot(fi0)=yal(fial);
    aad=ailus(fial,:);
    nad=nilus(fial,:);
    fival=find(aad>0);
    adu=ailuc(fi0,:);
    fival1=find(adu>0);
    aiadd=[ailuc(fi0,fival1(1:end-1))  aad(fival)];
    niadd=[niluc(fi0,fival1(1:end-1))  nad([fival fival(end)+1])];
    atot(fi0,1:length(aiadd))=aiadd;
    ftot(fi0,:)=fstus(fial,:);
    ntot(fi0,1:length(niadd))=niadd;
   end
   
   fiul=fial+1:length(ailus);
   fiulb=1:length(atot);
   sm=size(atot,2)-size(ailus,2)
   if sm>0
   zea=zeros(length(fiul),sm);
   zeb=[];
   else
   zeb=zeros(length(fiulb),-sm);
   zea=[];
   end
   atot=[[atot zeb]; [ailus(fiul,:) zea]];
   ntot=[[ntot zeb]; [nilus(fiul,:) zea]];   
   ftot=[ftot; fstus(fiul,:)];
   ytot=[ytot; yal(fiul)];
   
   figure, plot(atot,ytot,'c.')
pausak   
' da fare '
keyboard
ailua=atot;
nilua=ntot;
dilua=diff([0; ytot]);
fstua=ftot;

   fi=find(th<0);
   thb=abs(th(fi));
   nlb=[nl(1) nl(fi+1) nl(end)];
    
 if Np_ad~=0 
  [dilub,ailub,nilub,fstub]=lens_add(thb,h,r,N,nlb,ifps,Np_ad);
   fi0=find(ailub==0);
   ac_dis=ailub;
   ac_dis(fi0)=NaN;
   csub=cumsum(dilub);
 %  figure, plot(ac_dis,csu,'r.'), pausak

%   figure, plot(ac_dis,csub,'r.'), grid, pausak
 else
  dilub=[];
  ailub=[];
  nilub=[];
  fstub=[];
 end

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
   
   if Np_ad~=0
    fi=find(csub>h);
   else 
    fi=[];
   end 
%   figure, plot(ac_dis(fi,:),csua(fi),'g.'), grid, pausak
   figure, plot(ac_dis,csub,'r.',ac_dis(fi,:),csub(fi),'g.'), grid, pausak
   dilu1=[dilua; dilub(fi)];
   ailu1=[ailua; ailub(fi,:)];
   nilu1=[nilua; nilub(fi,:)];
   fstu1=[fstua; fstub(fi,:)];

   fi0=find(ailu1==0);
   a1_dis=ailu1;
   a1_dis(fi0)=NaN;
      csu1=cumsum(dilu1);
   figure, plot(a1_dis,csu1,'r.'),
   hold on
   fim=find(fstu1(:,2)>1);
   nrep=fstu1(fim([1 end]),2);
   
   col='cg';
   for k=1:length(nrep)
     fim=find(fstu1(:,2)==nrep(k));
     plot(a1_dis(fim,:),csu1(fim),[col(k),'.']),
   end
   title([' repetition ',num2str(nrep')])
   pausak

% per vedere lo specchi completo

 adisp=ailu1;
 ddisp=dilu1;
 ndisp=nilu1;
 fdisp=fstu1;
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
