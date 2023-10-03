function [fatqw,nref,aref,Lr,Lpl,nr,npl,xr,xpl,nmir,ar,shav,frp,lambda,Fi,any,...
          fmu,z,anys,uL,fiQ,fiC,Lpla,npla,Nspla,aral,...
          nv,dv,xm,radii,fst,perd,anyf,iauto,shavet,ipar,ifield,istfi,...
          icaplo,icsfit,icsfib,icsfi,istmet]=...
          parst_n(nomeF,lambda0,ifp,iany,ianys,nomeFe,Bi,ilo,icalcFi,imet);


if exist('icalFi')==0
 icalFi=1;
end
% read structure data

[nref,aref,nv,dv,xm,radii,fst,perd,anyf,iauto,dvw,xmw,fsw,dov,shavet,...
   ipar,ifield]=Layer_pa(lambda0,nomeF);

if length(find(ifield==-2))>1
 disp(' parst_pa.m: error: only one print section ')
 keyboard
end

disp('dopo layer_pa ')
keyboard
if length(find(ifield~=0))==0
 ifi=find(shavet~=0);
 ifield(ifi(1))=-2;
end

sip=size(ipar);
ip1=reshape(ipar(:,1,:),sip(1)*sip(3),1);
ip2=reshape(ipar(:,2,:),sip(1)*sip(3),1);
ip3=reshape(ipar(:,3,:),sip(1)*sip(3),1);
riga=[1:sip(1)]';
righe=riga*ones(1,sip(3));
ip4=reshape(righe,sip(1)*sip(3),1);
fi1=find(ip1~=0);
[dup,fi1p]=sort(ip1(fi1));
fi1=fi1(fi1p);

fid=find(diff([0; dup])~=0);
fiud=fi1(fid);
[ip11,fiudd]=sort(ip1(fiud));
fiu=fiud(fiudd);

npar=length(find(diff([0; dup])~=0));

 disp(' ')
 disp(' ')
 disp(' ')
 disp(' ')

screen('         # of parameters = ',npar)
 disp(' ')
 disp(' ')
 disp(' ')
type(1,:)='thickness    ';
type(2,:)='y-axis       ';
type(3,:)='x-axis       ';
type(4,:)='shape param  ';
type(5,:)='ref. index   ';
type(6,:)='array centers';
type(7,:)='shape        ';

types(1,:)='planar   ';
types(2,:)='rectangle';
types(3,:)='ellipse  ';
types(4,:)='rhombus  ';
types(5,:)='array    ';

%for k=1:length(fiu)
% screen('   --> parameter =',k)
% disp(['   --> layer =',num2str(ip4(fiu(k))),'   --> trans. section =',...
% num2str(ip3(fiu(k))), '  -->parameter type =',type(ip2(fiu(k)),:)]);
% disp(' ')
% pausak
%end

 disp(' ')
 disp(' ')
 disp(['      Parameters to be supplied in par_in ']);
 disp(' ')
 disp(['      Parameter,    Layer,   Trans. sect,   Par. type,     shape']);
 disp(' ')
 tab='          ';
 tab2='  ';
 tab1=' ';
 tab0='';
% fiu=fi1;

for k=1:length(fiu)
 tabi=tab1;
 if ip4(fiu(k))>=100
  tabi=tab0;
 elseif ip4(fiu(k))<10
  tabi=tab2;
 end
disp([tab num2str(ip1(fiu(k))) tab tabi  num2str(ip4(fiu(k))) tab ...
      num2str(ip3(fiu(k)))  tab type(ip2(fiu(k)),:)...
      tab2 types(shavet(ip4(fiu(k))),:)]);
end
 disp(' ')
 disp(' ')


%keyboard

 fiu=fi1;

%fis=find(shavet~=0);
%shd=shavet(fis);
%if length(find(diff(shd)~=0))>0
% istrumix=1;
%else
% istrumix=0;
%end
%istrumix=1;

reassign


xpl.b.last=xm(pup.b.last,1);
xpl.b.o=xm(pup.b.o,1);
xpl.b.m=xm(pup.b.m,1);
xpl.b.i=xm(pup.b.i,1);

xpl.t.o=xm(pup.t.o,1);
xpl.t.m=xm(pup.t.m,1);
xpl.t.i=xm(pup.t.i,1);



n1=(nv(:,1));
c1=xm(:,1);
xd=[zeros(size(dv)) dv];
fl=abs(fst(:,2));
fim=find(fl>1);
fw=length(fim);
ico=1;
piv=1;

dto=[];
nto=[];
ato=[];
cto=[];
puo=[];
fsto=[];
%'fw', keyboard
while fw>0
 pd=fim(ico);
 nl=fst(pd,1);
 nlay=abs(fst(pd,2));
 if piv-pd~=0
  pup0=piv:pd-1;
  if length(find(pup0==pua))==1
   fi0=find(pup0==pua);
   leqw=length(fsto)+fi0;
  end
  ay=anyf(pup0);
  nto=[nto; n1(pup0)];
  cto=[cto; c1(pup0)];
  fsto=[fsto; fst(pup0,2)];
  dto=[dto; dv(pup0)];
  ato=[ato; dv(pup0).*ay];
  puo=[puo; pup0'.*ay];
 end
 pupe=pd:pd+nl-1;
 ay=anyf(pupe);
 for inp=1:nlay
  nto=[nto; n1(pupe)];
  cto=[cto; c1(pupe)];
  fsto=[fsto; fst(pupe,2)];
  dto=[dto; dv(pupe)];
  if inp==1
   ato=[ato; dv(pupe).*ay];
   puo=[puo; pupe'.*ay];
  else
   ato=[ato; -dv(pupe).*ay];
   puo=[puo; pupe'*0];
  end
 end
 piv=pupe(length(pupe))+1;
 fw=fw-nl;
 ico=ico+nl;
end
if piv<length(dv)+1
 pup0=piv:length(dv);
 nto=[nto; n1(pup0)];
 cto=[cto; c1(pup0)];
 fsto=[fsto; fst(pup0,2)];
 dto=[dto; dv(pup0)];
 ay=anyf(pup0);
 ato=[ato; dv(pup0).*ay];
 puo=[puo; pup0'.*ay];
end

puac=puo;
puae=find(ato>0);


nd=[nto nto];
ccd=[cto cto];
xd=[zeros(size(dto)) dto];
xt=cumsum(reshape(xd',2*length(dto),1));
nt=real(reshape(nd',2*length(dto),1));
ct=(reshape(ccd',2*length(dto),1));

dJ=cumsum([0; dto]);
fJ=([0; fsto]);
nJ=([nto(1); nto]);
xJ=([cto(1); cto]);
lJ=length(dJ);
hz=[dJ(1):dJ(lJ)];
uFunc=hz*0;
uF0=hz*0;
perm=[];
xfit=[];
for il=1:lJ-1
 fiz=find(hz>=dJ(il) & hz<dJ(il+1));
 perm=[perm nJ(il+1)*ones(size(fiz))];
 xfit=[xfit xJ(il+1)*ones(size(fiz))];
 if fJ(il+1)==-1
  uFunc([fiz(1)-1 fiz])=1;
 end
 if fJ(il+1)==-1 & il==leqw
  uF0([fiz(1)-1 fiz])=1;
 end
end
xfit(length(hz))=xfit(length(hz)-1);
perm(length(hz))=perm(length(hz)-1);
uFunc(length(hz))=0;
uF0(length(hz))=0;
relPerm=conj(perm.^2);

 r411=-1.6e-6;  % in micron/V
 r412=-1.1e-6;
% r411=-2e-6;  % in micron/V
% r412=-2e-6;
 core=polyfit([0 1],[r411 r412],1);
 r41vet=polyval(core,xfit);
 n2r41=r41vet.*real(perm).^2;
 epsr41=2*r41vet.*real(perm).^4;
 fi0=find(xfit<0);
 n2r41(fi0)=0;
 epsr41(fi0)=0;

if iany>0
 wriwin
end

if icalcFi==1
 disp(' prima di risonanza ')
 keyboard

 [Ksi,lambda,Fi,uLong,uLong0]=eiglmio(lambda0,uFunc,uF0,relPerm);
 uL(1)=uLong0;
 uL(2)=uLong;
 disp(' fine Ksi '), keyboard
 fatqw=uLong/(uLong0*nmir.a);

else
 uL(1)=0;
 uL(2)=0;
 fatqw=0;
 lambda=0;
 Fi=0;
end

if ifp>=-1
 figure, plot(hz,real(perm),'w',hz,Fi*1.3,'r')
 figure, plot(xt,nt,'r.-',hz,real(perm),'w',hz,uFunc*4,hz,Fi)
 if ifp>1, pausak, end
end
%fatqw=1;
%Fi=1;
%lambda=lambda0;







if iany>0

 if ilo==0
  ico=0;
  for ifil=Bi
   ico=ico+1;
   fileName=[nomeFe,'_', num2str(ifil*10),'.dat'];
   disp([' loading file at V_bias = ',num2str(ifil)])
   [edu]=readsim(fileName);
   f(:,ico)=edu(:,2);
  end
  z=edu(:,1);
  parfe=' z f';
  eval(['save ' nomeFe parfe]);
 else
  eval(['load ' nomeFe]);
 end

 fmu=f*1e-4;
 finiz=find(nv(:,1)==1);
 dorigin=sum(dv(finiz));
 x=(hz'-dorigin)*1e-3;
 Fiany=Fi.*n2r41';
 Fianyeps=Fi.*epsr41';

 Fiwaeps=spline(x,Fianyeps,z);
 Fiwa=spline(x,Fiany,z);
 Fiw=spline(x,Fi,z);
 pewa=spline(x,imag(relPerm),z);

clear amedis ame amz amzd ameps

 dz=[diff(z); 0];
 smu=size(fmu);
 F_dz=((Fiwa.*dz)*ones(1,smu(2))).*fmu;
 F_dzeps=((Fiwaeps.*dz)*ones(1,smu(2))).*fmu;
 F_dis=((Fiw.*dz)*ones(1,smu(2))).*fmu;
 Fnorm=Fiw.*dz;
 Fnorms=Fiw'*dz;
 dmic=(dto)/1000;
 for k=1:length(puae)
  li=sum(dmic(1:puae(k)-1))-dto(1)/1000;
  lu=li+dmic(puae(k));
  fi=find(z>=li & z<lu);
  amedis(k,:)=sum(F_dis(fi,:))/sum(Fnorm(fi));
  ameps(k,:)=sum(F_dzeps(fi,:))/sum(Fnorm(fi));
  ame(k,:)=sum(F_dz(fi,:))/sum(Fnorm(fi));
  amz(k,:)=li;
  amzd(k,:)=dmic(puae(k));
 end

 iaa=1;
 amedis=amedis(:,iaa)';
 ad=[amedis amedis];
 xd=[amz' (amz+amzd)'];
 xat=(reshape(xd',2*length(amz),1));
 %xd=[zeros(size(amz')) amz'];
 %xat=cumsum(reshape(xd',2*length(amz),1));
 at=(reshape(ad',2*length(amz),1));

 xts=(xt-dorigin)/1000;
if ifp>=-1

 figure, plot(z,fmu,x,real(sqrt(relPerm)),'r',z,Fiw*5,'w',...
          xts,nt,'c',xts,nt,'m.',xat,at,'r.'),
 figure, plot(z,-fmu(:,iaa),'r',z,Fiw*50,'b',xts,nt.^2,'g'),
% figure, plot(z,fmu/100,x,sqrt(real(relPerm)),'y',z,Fiw*5,'w')
%          xts,nt*10,'c',xts,nt,'m.',xat,at,'r.'),
 if ifp>1, pausak, end
end

 anyv=zeros(length(dv),smu(2));
 puc=puac(find(puac>0));
 anyv(puc,:)=ameps;

% figure,
 puaf=[1:length(dto)];
 for k=1:length(puaf)
  li=sum(dmic(1:puaf(k)-1))-dto(1)/1000;
  lu=li+dmic(puaf(k));
  fi=find(z>=li & z<lu);
  if length(fi)==0
   Amedis(k)=0;
   Ame(k)=0;
  else
   sFn=sum(Fnorm(fi));
   if sFn~=0
    Amedis(k)=sum(F_dis(fi))/sFn;
    Ame(k)=sum(F_dz(fi))/sFn;
%    plot(z(fi),F_dis(fi),'.-',z(fi),Amedis(k)*ones(size(fi)),'w'), pausak
   else
    Amedis(k)=0;
    Ame(k)=0;
   end
  end
   Amz(k)=li;
   Amzd(k)=dmic(puaf(k));
 end

 Ad=[Amedis' Amedis'];
 Ade=[Ame' Ame'];
 Xd=[Amz' (Amz+Amzd)'-1e-6];
 Xat=(reshape(Xd',2*length(Amz),1));
 At=(reshape(Ad',2*length(Amz),1));
 Ate=(reshape(Ade',2*length(Amz),1));

 Atf=interp1(Xat,At,z,'linear');

 Atfe=interp1(Xat,Ate,z,'linear');


if ifp>=-1
 figure, plot(z,fmu,x,real(sqrt(relPerm)),'y',z,Fiw*5,'w',...
         xts,nt*10,'c',xts,nt,'m.',xat,at,'r.',Xat,At,'ro',z,Atf),
 if ifp>1, pausak, end
end


% return

 F_dz=Fiwa.*dz;
 ne=min(size(fmu));
 Fm=ones(ne,1)*F_dz';
 cSup=(Fm'.*fmu)/Fnorms;
 Sup=(F_dz'*fmu)/Fnorms;
 cobir=3e-1/lambda;
 birt=Sup*cobir

 fi0=find(pewa~=0);
 F_dz0=Fiwa.*dz;
 F_dz0(fi0)=0;
 Sup0=(F_dz0'*fmu)/Fnorms;
 bircav=Sup0*cobir
 fi0=find(pewa==0);
 F_dz0=Fiwa.*dz;
 F_dz0(fi0)=0;
 Sup0=(F_dz0'*fmu)/Fnorms;
 birmir=Sup0*cobir


% verifica

 F_dz=Atfe.*Fiw.*dz;
 ne=min(size(Atfe));
 Fm=ones(ne,1)*F_dz';
 Supv=sum(F_dz)/Fnorms;
 cSupv=(F_dz)/Fnorms;
 cobir=3e-1/lambda;
 birt_verifica=Supv*cobir

 if ifp>=-1
  figure, plot(z,cumsum(cSup)*cobir,z,cumsum(cSupv)*cobir,'-'), title(' fx-fy')
  if ifp>1, pausak, end
 end

 any.t=anyv(put,:);
 any.b=anyv(pub,:);
 any.a=anyv(pua,:);

else
 z=0;
 any=0;
 fmu=0;
end  %iany

if ianys>0
 p44=-0.072;
  anysv=-p44*ones(size(dv)).*real(nv(:,1)).^4;
  anys.t=anysv(put,:);
  anys.b=anysv(pub,:);
  anys.a=anysv(pua,:);
else
 anys=0;
end

%disp(' fine parstrut ')
%keyboard
