function [dilu,ailu,nilu,ipmir]=lens_sub(th,h,r,N,nl,updown,ifp,Rel,Nrel,Npair)
%function [dilu,ailu,nilu]=lens_sub(th,h,r,N,nl,updown,ifp,Rel,Nrel,Npair)
%' lens_sub in', keyboard
ipmir=0;
rin=r;
hin=h;
thin=th;
Nin=N;
if nargin<=9
 Npair=0;
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
 th=duth;
 dunr=[nl(1:nps+1) nl(end)];
 dusa0=dunr;
 dusa=dunr(2:end-1);
 clear dunr
 dunr=repmat(dusa,1,fix(Npair));
 if Npair-floor(Npair)>0
  dunr=[dunr dunr(1)];
 end
 dunr=[dusa0(1) dunr dusa0(end)];
 nl=dunr;
end

RMA=Rel;
if RMA==0
 RMA=20000;
end
%if updown==0
% r=-abs(r);
%end
if N==0
 dilu=[];
 ailu=[];
 nilu=[];
 return
end

%' sub in', keyboard
if N==1
 if h==0
  h=[];
 end
 dilu=[th h*1000]';
 ailu=zeros(size(dilu));
 nilu=nl(2:end-1)';
 if updown~=1
  dilu=flipud(dilu);
  ailu=flipud(ailu);
  nilu=flipud(nilu);
 end
%' sub in N=1', keyboard
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

hma0=min([h th]);
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
 xidu=sqrt(2*R*hv(n)-hv(n)^2);
 if xidu<RMA | k<=Nrel
  xi(ndu,k)=xidu;
 else
  xi(ndu,k)=RMA;
 end

% xi(ndu,k)=sqrt(2*R*hv(n)-hv(n)^2);
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
% if rm>RMA & k>Nrel
 if rm>RMA & k>=Nrel
  rm=RMA;
 end
 hmd=R-sqrt(R^2-rm^2);
 hm=hvero-hmd;
 if k>1000
% if k>000
 ' vedi ', keyboard
 end
 rim(n)=rm;
 him(n)=hm;

 A=Ai;
 riv=ri;
end
%if r>RMA & k>Nrel
if r>RMA & k>=Nrel
 r=RMA;
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
xi(1:length(him),k)=[fliplr(rim(2:end)) r]';


% ' vedi ', keyboard
for ke=k+1:length(hvet)
R1=Rvet(ke);
ci(:,ke)=R1*exp(j*tex);
hv1=h0-yi(:,k)+sum(th(1:ke-1));
for n=1:N+1
 yi(n,ke)=yi(n,k);
 xidu=sqrt(2*R1*hv1(n)-hv1(n)^2);
 if xidu<RMA | ke<=Nrel
  xi(n,ke)=xidu;
 else
  xi(n,ke)=RMA;
 end
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

%figure, plot(xl(:,k:end),yl(:,k:end),'--'), hold on, plot(xaud,yau,xi,yi,'r.'), pausak

end  %k

%return
%' qui ', pausak
xaus=xau;
%ic=1;
%for l=2:length(xau)
% if length(find(xau(l,:)~=0))==0
%  ic=ic+1;
%  xau(l,ic:end)=NaN;
% end
%end
 for kj=2:length(hvet)
  fim=find(xau(:,kj)==-1);
  if length(fim)>0
   xau(fim,kj)=xau(fim-1,kj);
  end
 end
 xl=real(ci);
 yl=imag(ci)-X0;

 if ifp==-10 | ifp>=0
 figure, plot(xl,yl,'--'), hold on, plot(xau,yau,xiu,yiu,'r.'),
 pausak
 end
%keyboard
%axis([0 r+th 0 h+th]),
dil=-diff(yau(1:2:end,end));
%dild=-[diff(yau(1:2:end,end)); 0];
%sai=size(xau);
%dild=repmat(dild,1,sai(2));
ail=(xau(2:2:end,:));
%' ail ', keyboard


fiz=find(dil~=0);


dilu=dil(fiz)*1000;
ailu=ail(fiz,:);
nilu=ones(length(ailu),1)*nl;

 keyboard
'lens_sub', keyboard

for kj=1:length(ailu)
 fi=find(isnan(ailu(kj,:))==1);
 if length(fi)>0
  fi1=find(isnan(ailu(kj,:))==0);
  ailud=ailu(kj,:);
  ailu(kj,:)=0;
  ailu(kj,1:length(fi1))=ailud(fi1);
  nilud=nilu(kj,:);
  nilu(kj,:)=0;
  fi1=[fi1 fi1(end)+1];
  nilu(kj,1:length(fi1))=nilud(fi1);
 end
end
% keyboard
%'lens_sub', keyboard

if rse==-1
 dilu=dilu(:,1);
 ailu=zeros(size(dilu(:,1)));
 nilu=nilu(:,1);
 inew=1;
 if inew==1
 dilu=[h0; th'];
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
%'lens_sub', keyboard
ailus=ailu;
for kk=1:length(dilu)
 du=ailus(kk,:);
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

if RMA==20000
 return
end
%          'lens_sub prima 1', keyboard

anu=zeros(size(ailu))*NaN;
nnu=zeros(size(nilu))*NaN;
for il=1:length(dilu)
% ddu=dilu(il,:);
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
 pun=[pua find(ndu==1)];
 anu(il,1:length(pua))=adu(pua);
 nnu(il,1:length(pun))=ndu(pun);
% pausak
end
adu=anu(:,1);
adus=adu;
dis=dilu;
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
fiM=find(adu==RMA);
isod=find(diff(fiM)==1);
fiso=fiM(isod);
adu(fiso)=-100;
z100=sum(dis(fiso))+dis(end);
fiR=find(adu~=-100);
anuu=anu(fiR,:);
nnuu=nnu(fiR,:);
dnuu=dis(fiR,:);
%dnuu(end)=z100;

anud=anuu;
fnN=find(isnan(anuu)==1);
anud(fnN)=-100;
lM=0;
for ka=1:length(anuu)
 fnN=length(find(anud(ka,:)~=-100));
 if fnN>lM
  lM=fnN;
 end
end
anuu=anuu(:,1:lM);
nnuu=nnuu(:,1:lM+1);
fi=find(isnan(anuu)==1);
anuu(fi)=0;
fi=find(isnan(nnuu)==1);
nnuu(fi)=0;
dilu=dnuu;
nilu=nnuu;
ailu=anuu;
if Npair<4
 fii=find(diff(flipud(nilu(:,1)))~=0);
% fii=find(diff((nilu(:,1)))~=0)+2;
 fiiu=fii([1 3]);
 %nilu(fiiu,:)
 %pausak
if updown==1
 ipmir=fiiu(1)+1:fiiu(end);
else
 ipmir=fiiu(1):fiiu(end)-1;
end
 dp=cumsum(flipud(dilu))/1000;
% dp=cumsum((dilu))/1000;
 od=ones(size(dp));

 if ifp==-10 | ifp>=0
  figure, plot(xl,yl,'--'), hold on,
  plot(xau,yau,xiu,yiu,'r.',od(ipmir)*0,dp(ipmir),'w.',...
  od(ipmir)*rin,dp(ipmir),'w.'),
  xlabel(' radial coordinate')
  ylabel(' longitudinal coordinate')
  grid
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

'fine lens_sub', keyboard
end
%'fine lens_sub', keyboard
