clear
close all
load spe
csa=sum(dilua);
fi=find(cumsum(dilub)>h);
figure, plot(cumsum(dilua),real(nilua(:,1)),'r.',csa+cumsum(dilub(fi)),real(nilub(fi,1)),'g.'),
hold on, plot(cumsum(dilua),real(fstua(:,2)),'r',csa+cumsum(dilub(fi)),real(fstub(fi,2)),'g'), pausak

pausak
%figure, plot(cumsum(dilua),real(nilua(:,1)),'r',csa+cumsum(dilub(fi)),real(nilub(fi,1)),'g'), pausak  

figure, plot(ailub(fi,:),csa+cumsum(dilub(fi)),'w.',ailua,cumsum(dilua),'.'), pausak

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
cs=cumsum(dilu);
figure, plot(cumsum(dilu),fstu(:,2),cs,nilu(:,1),'r.')
'qui cont', keyboard
isal=1;
if isal==0
dilu=[dilu(1); dilu];
nilu=[nilu(1,1)*ones(size(nilu(1,:))); nilu];
ailu=[ailu(1,:)*0; ailu];
fstu=[[0 1]; fstu];
ailu(end,:)=ailu(end,:)/2;

ailu=flipud(ailu);
nilu=flipud(nilu);
dilu=flipud(dilu)*1000;
fstu=flipud(fstu);
inew=1
if inew==1
if ifp<=-11
' modifica'
keyboard
end
end  % isal
ailu_s=ailu;
nilu_s=nilu;
%fi=find(ailu==r);
%ailu(fi)=0;
%%ailu(:,:)=0;
clear ailu
clear nilu
for k=1:length(ailu_s)
 fi=find(ailu_s(k,:)==r);
 if length(fi)>0
  if fi>1
   ailu(k,1:fi-1)=ailu_s(k,1:fi-1);
   nilu(k,1:fi)=nilu_s(k,1:fi);
  else
   ailu(k,1)=0;
   nilu(k,1)=nilu_s(k,1);
  end
 else
  ailu(k,:)=ailu_s(k,:);
  nilu(k,:)=nilu_s(k,:);
 end
end
sa=size(ailu);
k=sa(2);
 fi=find(ailu(:,k)==0);
 if length(fi)==sa(1);
  ailu=ailu(:,1:end-1);
  nilu=nilu(:,1:end-1);
 end

if ifp<=-11
' modifica'
keyboard
end
end


dilus=dilu;
ailus=ailu;
nilus=nilu;
fstus=fstu;
figure, plot(cumsum(dilu),nilu(:,1),'w.'), pausak
keyboard

fd=[1; fstu(:,2)];
pu=find(diff(fd)~=0);
fi=1
dilu=[];
ailu=[];
nilu=[];
pu=[pu; length(fstu)+1];
for k=1:length(pu)
 fu=pu(k)-1;
 fii=fi:fu;
[ fi fu fstu(fii(1),2)]
' fi ', pausak
figure, plot(cumsum(dilus(fii)),nilus(fii,1),'r.'), pausak
 fi=pu(k);

 dilu=[dilu; repmat(dilus(fii),fstu(fii(1),2),1)];
 ailu=[ailu; repmat(ailus(fii,:),fstu(fii(1),2),1)];
 nilu=[nilu; repmat(nilus(fii,:),fstu(fii(1),2),1)];
end 


dilud=[200/1000; dilu];
nilud=[1; nilu(:,1)];
dilud=[dilu];
nilud=[nilu(:,1)];
dd=repmat(dilud,1,2);
dd(:,1)=0;
dud=reshape(dd',length(dilud)*2,1);
nn=repmat(nilud,1,2);
nud=reshape(nn.',length(dilud)*2,1);
figure, plot(cumsum(dud),nud,'r.',cumsum(dud),nud,'w')
'in lensub_add ', keyboard
figure, plot(ailu,cumsum(dilu),'.'), pausak
