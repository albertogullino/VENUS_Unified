PUA=pua;

L1=[100 200]*1e3;
L3=[100 200]*1e3;

L1=[100:10:300]*1e3;
L2=[500:200:1500];
L3=[200]*1e3;

iEinc=1;
soloplot=1
iztutto=0;

for k1=1:length(L1)
for k2=1:length(L2)
for k3=1:length(L3)
 iauto(1,1)=0;
 di_ass=[L1(k1) L2(k2) L3(k3)]';
 dv(2:4)=di_ass;

  if iauto(1,1)==0
   fi=find(iauto(:,1)==1);
   iauto(fi,1)=0;
   iauto(1,1)=1;
   reassign
   Litn=Lr.t;
   nitn=nr.t.';
   aitn=ar.t.';

   fmlsp=abs(frp.t);
   fmlst=abs(frp.t);
   fi1=find(fmlsp(:,1)==0);
   fmlst(fi1,1)=1;
   fi2=find(fmlsp(:,2)==1);
   fmlst(fi1,2)=0;

   bi=aral.y.t';
   pai=aral.p.t';
   shai=shav.t';
   istfie=istfi.t;
   tyari=aral.ar.t';
   fmlsi=fmlst;

  end


  ifiez=2;
  freq=ze;
  [du,ifr]=min(abs(fre_camp-ze));
  ifr=ifr+1;
  ifr=1;
  kcav=kcav0*(1+freq);
  isaMo=1;

   tic
   acc0
   toc

%fie_new

  pua=PUA;
  [du,ilam]=min(abs(imag(Lamv(Nd1vu))));
  iN=Nd1vu(ilam);
  Lami=Lamv(iN);
  Ar=V(:,iN);
  A0=[Gas.*Ar; Ar];
  p=besm'*Ar;
  p2=abs(p).^2;
  ggg1=real(Lami)*c/(rr.*d*2);

  G00(iEinc)=ggg1/vgconv;
  E00(:,iEinc)=p2/max(p2);
  parinc(:,iEinc)=di_ass;


  if iztutto==1
   clear dmem Acoz
   ifiez=1;
   Lizi=5/1000;
   Liz0=50/1000;
   Liz0=2000/1000;

   re_new

   '[L1(k1) L2(k2) L3(k3)]'
   di_ass
  end  % iztutto
  iEinc=iEinc+1;
end
end
end
if iztutto==1
 figure, plot(xvero,E0)
else
 figure, plot(xvero,E00)
 pausak
 for k=1:length(parinc)
  plot(xvero,E00(:,k))
  text(max(xvero)/5,0.5,num2str(parinc(1,k)))
  text(2*max(xvero)/5,0.5,num2str(parinc(2,k)))
  text(3*max(xvero)/5,0.5,num2str(parinc(3,k)))
  pausak
 end

 for k=1:length(L2)
  lp=L2(k);
  fip=find(parinc(2,:)==lp);
  figure, plot(xvero,E00(:,fip))
  pausak
  figure
  for k1=fip
   plot(xvero,E00(:,k1)), pausak
  end
 end
end

  ifiez=0;
  return

%   tic
%   acc0
%   toc

fie_new

  ifiez=0;
%  sinte=size(Acoz);
  AF=Acoz(1:end/2,:);
  AB=Acoz(end/2+1:end,:);
  figure, plot(KK,abs(AF),KK,abs(Ar),'.-w')
%  figure, plot(KK,abs(AF),KK,abs(Afz),'.-w')
amn=max(AF(:,1));
amv=max(Afzf(:,1));
figure, plot(KK,abs(AF(:,[1 end])/amn)),
hold on, plot(KK,abs(Afzf(:,[1 end])/amv),'.')
figure, plot(KK,abs(AF(:,[1 end]))),
hold on, plot(KK,abs(Afzf(:,[1 end])),'.')
