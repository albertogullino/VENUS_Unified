%load pgam
%close all
%iLP=1;
Kp=0;
clear GGev GGmv Tev Tmv

nstv=[0:5:25];
nstv=[19];
nstv=Nspla.t.m;
freq=0;

%Lplit=[];
%nplit=[];
%Lplut=[];
%Lplut(1:3)=0;
%nplut=[];
%Lplt=[0.069; 0.060];
%nplt=[3.044; 3.5155];
%lambda=2*real(sum(Lplt.*nplt))
%pausak
ireal=1;
if ireal==1
 nplt=real(nplt);
 nplut=real(nplut);
 nplit=real(nplit);
end
dec=linspace(0,173,50);
dec=linspace(0,100,50);
Lplut0=Lplut;

for ist=1:length(nstv)
nstratiu=nstv(ist);

Lplut=Lplut0;
 for it=1:length(dec)
  dp=dec(it)*1e-3;
  if dp<=Lplut0(end)
   Lplut(end)=Lplut0(end)-dp;
  elseif dp>Lplut0(end) & dp<=sum(Lplut0(end-1:end))
   Lplut(end)=0;
   Lplut(end-1)=Lplut0(end-1)-dp+Lplut0(end);
  elseif dp>sum(Lplut0(end-1:end))
   Lplut(end)=0;
   Lplut(end-1)=0;
   Lplut(end-2)=Lplut0(end-2)-dp+sum(Lplut0(end-1:end));
  end
%  Lplut
%  dp
%  pausak


   [GGe1,GGm1,Te1,Tm1]=gaemms_h(Kp,freq,lambda,Lplit,nplit,Lplt,...
              nplt,nstratiu,rfu,rr,iLP1,Lplut,nplut,real(rqw));

   GGev(it,ist)=GGe1;
%  GGev(:,ist)=GGe1;
%  GGmv(:,ist)=GGm1;
%  Tev(:,ist)=Te1;
%  Tmv(:,ist)=Tm1;
 end
end
figure,
plot(dec,abs(GGev)),
