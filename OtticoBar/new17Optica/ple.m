%close all
clear
load qui   
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
    if fi0==402
     yyy
     yal(fial+[-1 0 1]), pausak
    end
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
   zea=zeros(length(fiul),size(atot,2)-size(ailus,2));
   atot=[atot; [ailus(fiul,:) zea]];
   zea=zeros(length(fiul),size(ntot,2)-size(nilus,2));
   ntot=[ntot; [nilus(fiul,:) zea]];   
   ftot=[ftot; fstus(fiul,:)];
   ytot=[ytot; yal(fiul)];
   
   figure, plot(atot,ytot,'c.')
