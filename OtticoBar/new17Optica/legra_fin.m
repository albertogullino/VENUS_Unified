clear 
close all
load grat
Np0=Npair;
rM=r;
%nl=nl.';
h=5e-2;
%thgv=[30 150];
thgv=[115.6 0];
%thgv=[.00001 80];
%thgv=[-th(end) 30 100];

%Np_ad=0;

if Np_ad==0
 Np_ad=-6.5;
 Npair=0;
 nl=[1 nl(4:6)];
 th=abs(th(3:4));
end


%Np0=-ceil(h*1000/sum(thgv))
%pausak
%if abs(Np0)>1
% Np0=Np0-1;
%end
if Npair~=0
  fi=find(th<=0);
else
  fi=find(th>=0);
end
fisla=fi;
    thb0=abs(th(fi));
    nlb0=(nl(fi+1));
    if abs(Np_ad)-floor(abs(Np_ad))>0
      nlb0=fliplr(nlb0);
    end  
    nlb=[ nlb0 -1  nl(end)];
%    nlb=[5 nlb0 -1  nl(end)];

fiz=find(thgv==0);
if length(fiz)==0
  thb=thgv;
else
  fizn=find(thgv~=0);
  thb=[thgv(fizn) 300];
  nlb=[nlb0 1  nl(end)];
end
    Np0=-2;
'ingresso', keyboard
   [dilub,ailub,nilub,fstub]=lens_add(thb,h,r,N,nlb,ifps,Np0);   

% prendo solo quello che interessa

% fidu=find(diff(fstub(:,2))~=0);
% fijou_in=fidu(1)+1;

'ingresso dopo', keyboard

Nsez=min([size(ailub,2) 3]);   
%fib=find(fstub(:,2)~=0);
fib=1:length(ailub);

%'Nzez', keyboard
ailub=ailub(fib,1:Nsez);
fi0=find(ailub==0);
ab_dis=ailub;
%ab_dis(fi0)=rM;
dilub=dilub(fib);
nilub=nilub(fib,1:Nsez+1);
fstub=fstub(fib,:);
      csub=cumsum(dilub);
%   figure, plot(ab_dis,csub,'r.'), pausak
   
R=(r^2+h^2)/(2*h);
R0=R-h;
y=R0+csub;
yy=repmat(y,1,size(ab_dis,2));
yy1=repmat(csub,1,size(ab_dis,2));
z=sqrt(yy.^2+ab_dis.^2)-R0;
rMax=h+sum(thgv/1000);
fiv=find(z<rMax*1.01);
five=find(z>rMax*1.01);
   figure, plot(ab_dis,csub,'r.',ab_dis(fiv),z(fiv),'w.'), pausak
   figure, plot(ab_dis(fiv),yy1(fiv),'c.'), pausak
   ap=ab_dis;
   ap(five)=0;
%   ap=ailub;
%   ap(five)=0;   
%  fijou=fijou_in:length(dilub);
  
   
    fi=fisla;   
    thb=abs(th(fi));
    nlb=[nl(1) nl(fi+1) nl(end)];
    Np0=sign(Np_ad)*(abs(Np_ad)-0.5);
%'ingresso', keyboard
   [diluc,ailuc,niluc,fstuc]=lens_add(thb,h,r,N,nlb,ifps,Np0);   
   fi0=find(ailuc==0);
   ac_dis=ailuc;
   ac_dis(fi0)=rM;
   csu=cumsum(diluc);
 %  figure, plot(ac_dis,csu,'r.'), pausak

   figure, plot(ac_dis,csu,'r.',ap,csub-h+csu(1)+csu(end),'c.'), grid
%   figure, plot(ac_dis,csu,'r.'),
   fire=find(fstuc(:,2)>1);
   hold on, plot(ac_dis(fire,:),csu(fire),'g.'), 
   title([' periodicita ',num2str(fstuc(fire(1),2))])
   pausak
   
%   fi0=find(ap(:,1)==0);
%   fiva=1:fi0(1)-1;
ymau=h+sum(thgv);

   fiva=find(csub<=ymau);
   ailus=ap(fiva,1:Nsez);
   nilus=nilub(fiva,1:Nsez+1);
   dilus=dilub(fiva);
   fstus=fstub(fiva,:);
   csubu=cumsum(dilus);
   yal=csubu-h+csu(end)-dilus(1);
   yal=csubu-h+csu(end);
   fia=find(ailuc(:,1)>0);
   pul=fia(end);
   yi=csu(pul)-h;
   yu=yi+h;
   ytot=csu; 
   [du,fi0]=min(abs(ytot-yi));
   yyy=ytot(fi0);
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

%    ytot(fi0)=yal(fial);
    aad=ailus(fial,:);
    nad=nilus(fial,:);
    fival=find(aad>0);
    adu=ailuc(fi0,:);
    fival1=find(adu>0);
    if length(fival1)>1
     aiadd=[ailuc(fi0,fival1(1:end-1))  aad(fival)];
     niadd=[niluc(fi0,fival1(1:end-1))  nad([fival fival(end)+1])];
    else
      if length(fival)>0
       aiadd=aad(fival);
       niadd=nad([fival fival(end)+1]);
      else
       if length(fival1)==1
        aiadd=ailuc(fi0,:);
        niadd=niluc(fi0,fival1(1:end-1));
       else 
        aiadd=[];
        niadd=[];
       end
      end
      
    end
    atot(fi0,1:length(aiadd))=aiadd;
    ftot(fi0,:)=fstus(fial,:);
    ntot(fi0,1:length(niadd))=niadd;
   end
   
   fiac=find(ailus(:,1)>0);
   fiul=fial+1:fiac(end);
   sa=size(atot,2);
   ss=size(ailus,2);
   if ss<sa
    zea=zeros(length(ailus),sa-ss);
    atot=[atot; [ailus(fiul,:) zea]];
    ntot=[ntot; [nilus(fiul,:) zea]];   
   else
    zea=zeros(length(atot),ss-sa);
    atot=[[atot zea]; ailus(fiul,:) ];
    ntot=[[ntot zea]; nilus(fiul,:) ];   
   end
   ftot=[ftot; fstus(fiul,:)];
   ytot=[ytot; yal(fiul)];
   
   figure, plot(atot,ytot,'c.')
pausak   

ailub=atot;
nilub=ntot;
dilub=diff([0; ytot]);
fstub=ftot;
 
if Npair~=0 
 
 fi=find(th>0);
 tha=th;
 nla=nl;
 [dilua,ailua,nilua,fstua]=lens_add(tha,h,r,N,nla,ifps,Npair);   
 
   fi0=find(ailua==0);
   ac_dis=ailua;
   ac_dis(fi0)=NaN;
   csua=cumsum(dilua);
 %  figure, plot(ac_dis,csu,'r.'), pausak

   figure, plot(ac_dis,csua,'r.'), grid, pausak


   sa=size(ailua,2);
   sb=size(ailub,2);
   if sa~=sb
    if sa>sb
     le=length(dilub);
     sm=sa-sb;     
     ze=zeros(le,sm);
     nilub=[nilub ze];
     ailub=[ailub ze];
    else
     le=length(dilua);
     sm=sb-sa;
     ze=zeros(le,sm);
     nilua=[nilua ze];
     ailua=[ailua ze];
    end
   end
   
   if Np_ad~=0
    fidu=find(diff(fstub(:,2))~=0);
    fi=fidu(1)+1:length(fstub);
   else 
    fi=[];
   end 
   figure, plot(ac_dis,csua,'r.'), grid, pausak
   dilu1=[dilua; dilub(fi)];
   ailu1=[ailua; ailub(fi,:)];
   nilu1=[nilua; nilub(fi,:)];
   fstu1=[fstua; fstub(fi,:)];
else
   dilu1=[dilub];
   ailu1=[ailub];
   nilu1=[nilub];
   fstu1=[fstub];
end

fi=find(dilu1>1e-5);
 adisp1=ailu1(fi,:);
 ddisp1=dilu1(fi);
 ndisp1=nilu1(fi,:);
 fdisp1=fstu1(fi,:);
fi=find(ndisp1(:,1)~=0);
 adisp=adisp1(fi,:);
 ddisp=ddisp1(fi);
 ndisp=ndisp1(fi,:);
 fdisp=fdisp1(fi,:); 
 
 dilu1=ddisp;
 ailu1=adisp;
 nilu1=ndisp;
 fstu1=fdisp;

   fi0=find(ailu1==0);
   a1_dis=ailu1;
   a1_dis(fi0)=NaN;
      csu1=cumsum(dilu1);
   figure, plot(a1_dis,csu1,'r.'),
   hold on
   fim=find(fstu1(:,2)>1);
   df=diff(fim);
   if length(find(df>1))>0
    nrep=fstu1(fim([1 end]),2);
   else
    nrep=fstu1(fim(1),2);
   end    
   
   col='cg';
   for k=1:length(nrep)
     fim=find(fstu1(:,2)==nrep(k));
     plot(a1_dis(fim,:),csu1(fim),[col(k),'.']),
   end
   title([' repetition ',num2str(nrep')])
   pausak

% per vedere lo specchi completo

% adisp=ailu1;
% ddisp=dilu1;
% ndisp=nilu1;
% fdisp=fstu1;
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

  
    figure,  plot(dd0,nn0,csu1(fc-1),n0(fc-1),'.'), pausak

figure, plot(cumsum(dilu1),nilu1,'.')