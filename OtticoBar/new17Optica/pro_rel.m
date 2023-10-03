
  Lp=flipud(Litn);
  np=flipud(nitn(1,:).');
  flp=flipud(fmlst);
  %fi1=find(flp(:,1)==1);
  du=diff([flp(:,1)]);
  fii0=find(du~=0);
  fii=[1:fii0(1)];
  Lpi=Lp(fii);
  npi=np(fii);
  fim=find(flp(:,1)>1);
  Lpm=Lp(fim);
  npm=np(fim);
  npair=flp(fim(1),2);

  fiiu=[fii0(2)+1:length(Lp)];
  Lpu=Lp(fiiu(1));
  npu=np(fiiu(1));
  Lpu0=sum(Lpu);
  npu1=[nv(4,1)];
  npu2=[nv(2,2)];
  npu3=[nv(3,2)];
  [ndu1,ndu2]=nti(lambda*1e-6);
  npu3p=ndu1-j*ndu2;
  Lpu2=[dv(2)]*1e-3;
  Lpu3=[dv(3)]*1e-3/3;
  disp(' pro_rel ')
  keyboard

  spac=linspace(0,.2,200);
%  kv=[0 .1];
clear Gape Game
clear GapeM GameM
  kv=0;
%  for ks=1:length(spac)
%   spaci=spac(ks);
%   Lpu=[spaci Lp(end)];
%   npu=npu1;
%   [Gae,Gam]=gaemms(kv,0,lambda,Lpi,npi,Lpm,npm,npair,rfu,rr,iLP,Lpu,npu);
%   GapeM(ks,:)=(Gae).';
%   GapmM(ks,:)=(Gam).';
%  end

  kv=0;
  for ks=1:length(spac)
   spaci=spac(ks);
   Lpu=[spaci];
   npu=npu1;
   [Gae,Gam]=gaemms(kv,0,lambda,Lpi,npi,Lpm,npm,npair,rfu,rr,iLP,Lpu,npu);
   Gape(ks,1)=(Gae).';
   Gapm(ks,:)=(Gam).';
  end
  kv=0;
  for ks=1:length(spac)
   spaci=spac(ks);
   Lpu=[spaci Lpu3 Lpu3 Lpu3 Lpu2];
   npu=[npu1 npu3 npu3p npu3 npu2];
   [Gae,Gam]=gaemms(kv,0,lambda,Lpi,npi,Lpm,npm,npair,rfu,rr,iLP,Lpu,npu);
   Gape(ks,2)=(Gae).';
   Gapm(ks,:)=(Gam).';
  end
%  figure, plot(spac,abs(Gape)), hold on, plot(spac,abs(Gapm),'--'),
%  figure, plot(spac,abs(Gape),spac,abs(GapeM),'--'),
%  title(' yell no met, magenta: met')
  figure, plot(spac*1000,abs(Gape)),
% axis([0 200 .97 1])
%  title(' yell no met, magenta: met')
