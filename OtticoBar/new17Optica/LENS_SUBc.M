function [dilu,ailu,nilu,fstu]=lens_sub(th,h,r,N,nl,updown,ifp,Rel,Hre,Npair,Rflat,Rm_ring)

'in lensub ', keyboard
Rmax=r;
if nargin<=11
 Nrel=Hre;
else
% if Rm_ring*Hre>0
 if Rm_ring>0 & (N>1 | N<0)
  [dilu,ailu,nilu,fstu]=lens_rel(th,h,r,N,nl,ifp,Hre,Npair,Rel,Rm_ring);
%  [dilu,ailu,nilu,fstu]=lens_sav(th,h,r,N,nl,ifp,Hre,Npair,Rel,Rm_ring);
  return
 else
  Nrel=0;
 end
end

if Nrel>0
 Nrels=(2*Npair-Nrel)+1;
elseif Nrel<0
 Nrels=Nrel-1;
else
 Nrels=0;
end
%'in lens_sub dope Nrels', keyboard

%ifp=-11
ipmir=0;
nin=nl;
rin=r;
hin=h;
ha=h;
thin=th;
Nin=N;
Npa_in=Npair;
duth=th;
' lens_sub in 0', keyboard
if nl(1)==nl(end)
 NPadd=1;
else
 NPadd=2;
end
% NPadd=1;

  if Npair>0
   no_pair=1;
   no_rep=0;
  else
   per0=sum(th(1:2))/1000;
%   Npe=floor(2*ha/per0)+NPadd;
%   Npe=floor(2*ha/per0)+1 + abs(Nrel/2);
   Npe=ceil(ha/per0)+1 + abs(Nrel/2);
   if Npair-floor(Npair)>0
    Npadu=Npe+.5;
   else
    Npadu=Npe;
   end
   no_pair=floor(abs(Npair))-Npe+1;
%   no_pair=floor(abs(Npair))-Npe;
   Npair=Npadu;
  end
%' lens_sub in', keyboard
if nargin<=9
 Npair=0;
end
if N==1
 Npair=abs(Npa_in);
 Npair=1;
 Npair=abs(Npa_in)-fix(abs(Npa_in))+1;
end
if Npair>0
 nps=2;
 duth=th(1:nps);
 dusa=duth;
 clear duth
 duth=repmat(dusa,1,fix(Npair));
 if Npair-floor(Npair)>0
  duth=[duth duth(1)];
 end
 if length(thin)>nps
  duth=[duth thin(nps+1)];
 end
 th=duth;
 dunr=[nl(1:nps+1) nl(end)];
 dusa0=dunr;
 dusa=dunr(2:end-1);
 clear dunr
 dunr=repmat(dusa,1,fix(Npair));
 if Npair-floor(Npair)>0
  dunr=[dunr dunr(1)];
 end
 if length(thin)>nps
  dunr=[dunr nin(nps+2)];
 end
 dunr=[dusa0(1) dunr dusa0(end)];
 nl=dunr;

end


%' dopo', keyboard
RMA1=Rel+Rflat;
RMA=Rel;
if RMA==0
 RMA=20000;
end
%if updown==0
% r=-abs(r);
%end
%' sub in', keyboard
if N==0
 dilu=[];
 ailu=[];
 nilu=[];
 fstu=[];
 return
end

%' sub in', keyboard
if N==1
 if h==0
  h=[];
  nilu=nl(2:end-1)';
  dilu=th(2:end-1)';
 else
%  nilu=nl(1:end-1)';
  nilu=nl(1:end-1)';
 end
 dilu=[h*1000 th]';
 ailu=zeros(size(dilu));
 if updown~=0
  dilu=flipud(dilu);
  ailu=flipud(ailu);
  nilu=flipud(nilu);
 end
 ipmir=1+[1:nps];
% ' ipmir', keyboard
  if Npa_in>0
   no_pair=1;
   no_rep=0;
   pu_fst=[ones(size(dilu))*no_rep ones(size(dilu))*no_pair];
  else
   no_rep=length(ipmir);
   pu_fst=[ones(size(dilu))*0 ones(size(dilu))*1];
   pu_fst(ipmir,1)=no_rep;
   pu_fst(ipmir,2)=no_pair;
  end
  fstu=pu_fst;
  if ifp==-10
%   ' sub in N=1', keyboard
  end
 return
end

%' sub in', keyboard
rse=sign(r);
r=abs(r);
th=th/1000;
if ~exist('ifp')
 ifp=-10;
end

%h=R-sqrt(R^2-r^2);
h0=h;
R=(r^2+h^2)/(2*h);
X=(R-h);
X0=X;
te=atan(r/X);
tex=pi/2-linspace(0,te,100)';
%figure, polar(tex,R*ones(size(tex)))
Nx=100;
x=linspace(0,r,Nx)';
hveto=[h h+cumsum(th)];
hvetou=[0 hveto];
hvet=[h th];

%keyboard
Rvet=[R R+cumsum(th)];
Xvet=Rvet-hvet;
rvet=sqrt(2*hvet.*Rvet-hvet.^2);
tevet=atan(rvet./Xvet);
%keyboard
%
% rdisc=r/N;
% rv=[rdisc:rdisc:r];
% hdisc=R-sqrt(R^2-rv.^2);
%
%
%
%
%
%
%
%


xau=[];
yau=[];
xiu=[];
yiu=[];

if Rflat==0
 hma0=max([h th]);
else
 hma0=min([h th]);
end
dh0=hma0/N;

hvet=hvet;
hvet1=[0 hveto];



ndu=0;

for k=1:length(hvet)
%for k=1:2
%for k=2:2
%for k=3:3
%for k=2:3
%for k=1:1
ndu=0;
clear xa ya rim him
clear xi yi
hma=hvet(k);
if hma==hma0
 dh=hma/N;
else
 dh=diff(hvet1(k:k+1));
 N=ceil(dh/dh0);
 dh=dh/N;
end
hv=[hvet1(k):dh:hvet1(k+1)]-hvet1(k);
hv1=[hvet1(k):dh:hvet1(k+1)];
%pausak


if k>1
% hv=hv(2:end);
% N=N-1;
end
%keyboard
R=Rvet(k);
h=hvet(k);
hvero=hveto(k);
r=rvet(k);
ci(:,k)=R*exp(j*tex);
for n=1:N+1
 ndu=ndu+1;
 ndu=n;
 yi(ndu,k)=hvero-hv(n);
% [ndu k]
% pausak
 xidu=sqrt(2*R*hv(n)-hv(n)^2);
 if Nrels>0
  if xidu<RMA | k<=Nrels
   xi(ndu,k)=xidu;
  else
   xi(ndu,k)=RMA;
  end
 elseif Nrels<0
  if xidu<RMA | k>=abs(Nrels)
   xi(ndu,k)=xidu;
  else
   xi(ndu,k)=RMA;
  end
 elseif Nrels==0
  if xidu<RMA
   xi(ndu,k)=xidu;
  else
   xi(ndu,k)=RMA;
  end
 end
end

riv=r;
X=(R-h);
te=atan(r/X);
l=te*R;
A=R*l-r*X;
hvero=hveto(k);
for n=2:N+1
 hi=hv1(n);
 hid=hvero-hi;
 Xi=(X0+hi);
 ri=sqrt(2*R*hid-hid^2);
 tei=atan(ri/Xi);
 li=tei*R;
 Ai=R*li-ri*Xi;
 dA=A-Ai;
 rm=dA/dh/2;
 isa=0;
 if isa==0
 rms=rm;
 if Nrels>0
  if rm>RMA & k>=Nrels
   rm=RMA;
  end
 elseif Nrels<0
  if rm>RMA & k<abs(Nrels)
   rm=RMA;
  end
 elseif Nrels==0
  if rm>RMA
   rm=RMA;
  elseif max(rm)==0
   'rm=0', keyboard
  end
 end
 end %isa
 hmd=R-sqrt(R^2-rm^2);
 hm=hvero-hmd;
 hmds=R-sqrt(R^2-rms^2);
 hms=hvero-hmds;
 if k>1000
% if k>000
% ' vedi ', keyboard
 end
 rim(n)=rm;
 if rm~=rms
  him(n)=hms;
 else
  him(n)=hm;
 end
% [hm hms]
% 'hm', pausak

 A=Ai;
 riv=ri;
end

 if Nrels>0
  if r>RMA & k>=Nrels
   r=RMA;
  end
 elseif Nrels<0
  if r>RMA & k<=abs(Nrels)
   r=RMA;
  end
 elseif Nrels==0
  if r>RMA
   r=RMA;
  end
 end
rimu=fliplr(rim);
rimu=[0 rimu(1:end-1) r];
rd=[rimu(1:end-1); rimu(2:end)];
hd=hvero-[hv(1:end); hv(1:end)];
rp=reshape(rd,1,2*length(rd));
hp=reshape(hd,1,2*length(rd));

xa(:,k)=rp';
ya(:,k)=hp';


hma=h;


yi(1:length(him),k)=[fliplr(him(2:end)) hvetou(k)]';
%' yi ', keyboard
xi(1:length(him),k)=[fliplr(rim(2:end)) r]';

for ke=k+1:length(hvet)
R1=Rvet(ke);
ci(:,ke)=R1*exp(j*tex);
hv1=h0-yi(:,k)+sum(th(1:ke-1));
for n=1:N+1
 yi(n,ke)=yi(n,k);
 xidu=sqrt(2*R1*hv1(n)-hv1(n)^2);
 if Nrels>0
  if xidu<RMA | ke<=Nrels
   xi(n,ke)=xidu;
  else
   xi(n,ke)=RMA;
  end
 elseif Nrels<0
  if xidu<RMA | ke>=abs(Nrels)
   xi(n,ke)=xidu;
  else
   xi(n,ke)=RMA;
  end
 elseif Nrels==0
  if xidu<RMA
   xi(n,ke)=xidu;
  else
   xi(n,ke)=RMA;
  end
 end
% [n k]
 ylay(n,k)=k;

end

rimu=[0 xi(:,ke)'];
rd=[rimu(1:end-1); rimu(2:end)];
hd=hvero-[hv(1:end); hv(1:end)];
rp=reshape(rd,1,2*length(hv));
hp=reshape(hd,1,2*length(hv));

xa(:,ke)=rp';
ya(:,ke)=hp';
%' ke', keyboard

end  %ke
 xa(:,1:k-1)=NaN;
 xa(1,k+1:end)=NaN;
 xa(1,k+1:end)=-1;
%[xa ya]
xau=[xa; xau];
yau=[ya; yau];
xid=xi;
if k>1
 xid(end,:)=NaN;
end
xiu=[xid; xiu];
yiu=[yi; yiu];

xl=real(ci);
yl=imag(ci)-X0;
  fim=find(xau==-1);
  xaud=xau;
  xaud(fim)=NaN;


end  %k

%' quo sub'
%keyboard
xaus=xau;
 for kj=2:length(hvet)
  fim=find(xau(:,kj)==-1);
  if length(fim)>0
   xau(fim,kj)=xau(fim-1,kj);
  end
 end
 xl=real(ci);
 yl=imag(ci)-X0;

 if ifp<=-11 | ifp>=0
 figure, plot(xl,yl,'--'), hold on, plot(xau,yau,xiu,yiu,'r.'),
 pausak
 end
% dil=-diff(yau(1:2:end,end));

%clear xauu yauu
%for kk=1:length(xau)
% fi=find(isnan(xau(kk,:))==0 & xau(kk,:)~=0 );
% xauu(kk,1:length(fi))=xau(kk,fi);
% yauu(kk,1:length(fi))=yau(kk,fi);
%end
xauu=xau;
yauu=yau;


% dil=-diff(yau(1:2:end,1));
% ail=(xau(2:2:end,:));

 dil=-diff(yauu(1:2:end,end));
 ail=(xauu(2:2:end,:));
%' ail ', keyboard



fiz=find(dil>1e-7);


dilu=dil(fiz)*1000;
ailu=ail(fiz,:);
nilu=ones(length(ailu),1)*nl;
%' nilu ', keyboard


for kj=1:length(ailu)
 fi=find(isnan(ailu(kj,:))==1);
 if length(fi)>0
% kj, 'nilui i ', pausak
  fi1=find(isnan(ailu(kj,:))==0);
  ailud=ailu(kj,:);
  ailu(kj,:)=0;
  ailu(kj,1:length(fi1))=ailud(fi1);
  nilud=nilu(kj,:);
  nilu(kj,:)=0;
  fi1=[fi1 fi1(end)+1];
  nilu(kj,1:length(fi1))=nilud(fi1);
%  real(nilu(kj,:)),
%  ailu(kj,:)
%  'nilui u ', pausak
  end
end
% keyboard
%'lens_sub nilui', keyboard

if rse==-1
 dilu=dilu(:,1);
 ailu=zeros(size(dilu(:,1)));
 nilu=nilu(:,1);
 inew=1;
 if inew==1
 dilu=[h0; th'];
 if length(find(dilu<0))>0
  '<0 ', keyboard
 end
 ailu=zeros(size(dilu(:,1)));
 nilu=nl(1:end-1)';
 dilu=flipud(dilu)*1000;
 ailu=flipud(ailu);
 nilu=flipud(nilu);
% keyboard
 end

end

if updown~=1
 dilu=flipud(dilu);
 ailu=flipud(ailu);
 nilu=flipud(nilu);
end
%'lens_sub Rflat', keyboard

if Rflat>0
 fi=find(ailu~=0);
 ailu(fi)=ailu(fi)+Rflat;
end
ailus=ailu;
for kk=1:length(dilu)
 du=ailus(kk,:);
% fiM=find(du==RMA1);
 fiM=find(du==RMA);
 if length(fiM)>0
  fiR=fiM(1);
 else
  fiR=0;
 end
 if fiR>1
  du=du(fiR-1:fiR);
  if diff(du)<0
   ailu(kk,fiR:end)=0;
  end
 end
end
' dopo RMA ', keyboard
if RMA==20000
 return
end
%          'lens_sub prima 1', keyboard
isa=0
if isa==0
anu=zeros(size(ailu))*NaN;
nnu=zeros(size(nilu))*NaN;
for il=1:length(dilu)
 adu=ailu(il,:);
 ndu=nilu(il,:);
 fi=find(adu==RMA);
 if length(fi)==0
  fi=find(adu==0)-1;
  if length(fi)==0
   fi=length(adu);
  end
 end
 pua=1:fi(1);
% pun=[pua find(ndu==1)];
 pun=[pua find(ndu==nin(end))];
 anu(il,1:length(pua))=adu(pua);
 nnu(il,1:length(pun))=ndu(pun);
% pausak
end
end

adu=anu(:,1);
adus=adu;
dis=dilu;
%' dis', keyboard
su=0;
for ka=1:length(adus)
 adua=adus(ka);
 if adua==RMA
  su=su+dilu(ka);
 else
  if su>0
   dis(ka-1)=su;
   su=0;
  end
 end
end

%fiM=find(adu==RMA1);
%isod=find(diff(fiM)==1);
%fiso=fiM(isod);
%adu(fiso)=-100;
%z100=sum(dis(fiso))+dis(end);
%fiR=find(adu~=-100);
fiR=1:length(anu);
anuu=anu(fiR,:);
nnuu=nnu(fiR,:);
dnuu=dis(fiR,:);
%dnuu(end)=z100;
%' controlo', keyboard

anud=anuu;
fnN=find(isnan(anuu)==1);
anud(fnN)=-100;
lM=0;
sa=size(anuu);
for ka=1:length(anuu)
 fnN=length(find(anud(ka,:)~=-100 & anud(ka,:)~=RMA1));
 if fnN>lM & fnN<sa(2)
  lM=fnN+1;
 end
end
anuu=anuu(:,1:lM);
nnuu=nnuu(:,1:lM+1);

anuuu=[];
nnuuu=[];
dnuuu=[];
kk=1;
while kk<=length(anuu)
%kk
 if anuu(kk,1)==RMA1
%  ' passo '
  fi=find(anuu(kk:end,1)==RMA1);
  dfi=find(diff(fi)>1);
  if length(dfi)>0
   puu=dfi(1)-1;
  else
   puu=length(fi)-1;
  end
  pu=kk+[0:puu];
  lpu=length(pu);
%  keyboard
  anuuu=[anuuu; anuu(kk,:)];
  nnuuu=[nnuuu; nnuu(kk,:)];
  dnuuu=[dnuuu; lpu*dnuu(kk,:)];
  kk=kk+lpu;
 else
%  ' passo 0'
  anuuu=[anuuu; anuu(kk,:)];
  nnuuu=[nnuuu; nnuu(kk,:)];
  dnuuu=[dnuuu; dnuu(kk,:)];
  kk=kk+1;
 end
% pausak
end
anuus=anuu;
anuu=anuuu;
dnuus=dnuu;
dnuu=dnuuu;
nnuus=nnuu;
nnuu=nnuuu;
%' cont qwuo', keyboard

fi=find(isnan(anuu)==1);
anuu(fi)=0;
fi=find(isnan(nnuu)==1);
nnuu(fi)=0;
dilu=dnuu;
nilu=nnuu;
ailu=anuu;
if Npair<20
 if updown==1
  fii=find(diff(flipud(nilu(:,1)))~=0);
%  fiid=find(diff((nilu(:,1)))~=0);
  fiid=find(diff(flipud(nilu(:,1)))~=0);
%  dp=cumsum(flipud(dilu))/1000;
  dp=cumsum((dilu))/1000;
 else
  fii=find(diff((nilu(:,1)))~=0);
%  fiid=find(diff(flipud(nilu(:,1)))~=0);
  fiid=find(diff((nilu(:,1)))~=0);
  dp=cumsum((dilu))/1000;
%  dp=cumsum(flipud(dilu))/1000;
 end
% 'fii', keyboard
  fii=find(diff((nilu(:,1)))~=0);
%  fiid=find(diff(flipud(nilu(:,1)))~=0);

 %nilu(fiiu,:)
% 'fiid'
% pausak
   fiiud=fiid([1 3]);
  if updown==0
   fiiu=fii([1 3]);
  else
   fiiu=fii([end-2 end]);
  end
   ipmir=fiiu(1)+1:fiiu(end);
%  pausak
  ipmird=fiiud(1)+1:fiiud(end);
 od=ones(size(dp));

 if ifp<=-11 | ifp>=0
  figure, plot(xl,yl,'--'), hold on,
  plot(xau,yau,xiu,yiu,'r.',od(ipmird)*0,dp(ipmird),'w.',...
  od(ipmird)*rin,dp(ipmird),'w.'),
  xlabel(' radial coordinate')
  ylabel(' longitudinal coordinate')
  grid
  pausak
 end
 dperdu=dp(ipmir);
 dper=diff(dperdu([1 end]));
 dperb=dp(ipmir(1));
% ipmirb=fiiu(end)+1:length(dilu);
% dperdub=dp(ipmirb);
% dperb=diff(dperdub([1 end]));
else
 ipmir=0;
end


if ifp~=-4
%figure, plot(anuu,cumsum(dnuu),'.'), hold on, plot(xau,yau*1000,'o')
%'fine lens_sub', keyboard
end
%ailu=ailu+Rflat;
fi=find(ailu==Rel);
%ailu(fi)=0;
%if Rflat>0
% fi=find(ailu~=0);
% ailu(fi)=ailu(fi)+Rflat;
%end
%' Rel ', keyboard

  if Npa_in>0
   no_pair=1;
   no_rep=0;
   pu_fst=[ones(size(dilu))*no_rep ones(size(dilu))*no_pair];
  else
   no_rep=length(ipmir);
   pu_fst=[ones(size(dilu))*0 ones(size(dilu))*1];
   pu_fst(ipmir,1)=no_rep;
   pu_fst(ipmir,2)=no_pair;
  end
fstu=pu_fst;

%if Rflat>0
%if Rflat<-1
if Rflat>-1
 ru=Rflat+r;
 for k=1:length(dilu)
  ad=ailu(k,:);
  fia=find(ad==ru);
  if length(fia)>=1
   fie=fia(1);
   ailu(k,fie:end)=0;
   nilu(k,fie:end)=nilu(k,fie);
  end
 end
 sa=size(ailu);
 for k=1:sa(2)
  fi=find(ailu(:,k)==0);
  if length(fi)==sa(1)
   ize(k)=1;
  else
   ize(k)=0;
  end
 end
 fiza=find(ize==0);
% fizn=[fiza sa(2)+1];
 fizn=[fiza fiza(end)+1];
 ailus=ailu;
 nilus=nilu;
 ailu=ailu(:,fiza);
 nilu=nilu(:,fizn);
end


itro=find(diff(fstu(:,2))~=0);
itro=[itro; length(fstu)];
ki=1;
xd=[];
yd=[];
dsu=cumsum(dilu)/1000;
sup=0;
for k=1:length(itro)
 ku=itro(k);
 pu=ki:ku;
 ydu=dsu(pu);
 xdu=ailu(pu,:);
 RP=fstu(ki,2);
 if RP==0
  RP=1;
 end
 for kk=1:RP
  if kk>1
   sup=sup+per0;
  end
  xd=[xd; xdu];
  yd=[yd; ydu+sup];
 end
% if kk<RP
%  sup=sup+per01;
%  xd=[xd; xdu];
%  yd=[yd; ydu+sup];
% end
 ki=ku+1;
end
if ifp<=-11
 figure
 plot(xd,yd,'.'), grid
 pausak
 figure, plot(fstu(:,2),'.')
 pausak

end
%fi=find(ailu==0);
%ailu(fi)=Rmax;


if ifp<=-10
'fine lens_sub',
%pausak
keyboard
end
