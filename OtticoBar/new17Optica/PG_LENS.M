load le_and
%load pg_old
%close all
%iLP=1;
Kp=[0:.05:.25];
Kpf=[0:.005:.3];

Kp=[0];
Kpf=[0];
clear GGev GGmv Tev Tmv GGek GGmk

nstv=[0:5:25];
nstv=Nspla.t.m;

%nstv=[23];

%Lplit=[];
%nplit=[];
%Lplut=[];
%Lplut(1:3)=0;
%nplut=[];
%Lplt=[0.069; 0.060];
%nplt=[3.044; 3.5155];
%lambda=2*real(sum(Lplt.*nplt))
%pausak

fi=find(Lplut>0 & nplut(:,1)>2);
%fi=[fi; fi(end)+1];
Lplut=Lplut(fi);
nplut=nplut(fi,:);

%keyboard
Lplut(end)=102/1000;

%Lplut(end-1)=38/1000;
%Lplut(end)=90/1000;
%
%Lplut(end-1)=0/1000;
%Lplut(end-1)=60/1000;
%Lplut(end)=113/1000;

% old
%Lplut(end-1)=0/1000;
%Lplut(end)=120/1000;


%keyboard

ireal=1;
if ireal==1
 nplt=real(nplt);
 nplut=real(nplut);
 nplit=real(nplit);
end

nse=1;
ico=0;

ma_dec=sum(Lplut(end-nse:end))*1000;
dec=linspace(0,ma_dec,50);
dec=linspace(0,ma_dec,100);
dec=linspace(0,80,100);
%dec=linspace(0,80,81);
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
              nplt,nstratiu,rfu,rr,0,Lplut,nplut,real(rqw));

%   GGev(it,ist)=GGe1;
   GGev(it,:)=GGe1.';
   GGmv(it,:)=GGm1.';

  if floor(it/10)-ico~=0 & length(Kpf)>1
   ico=ico+1;
   [GGe1,GGm1,Te1,Tm1]=gaemms_h(Kpf,freq,lambda,Lplit,nplit,Lplt,...
              nplt,nstratiu,rfu,rr,0,Lplut,nplut,real(rqw));

   GGek(ico,:)=GGe1.';
   GGmk(ico,:)=GGm1.';
  end

%  Tev(:,ist)=Te1;
%  Tmv(:,ist)=Tm1;
 end
end


figure
T=1-abs(abs(GGev)).^2;
plot(dec,T),

return
pausak
figure
plot(dec,abs(GGev)),

return
figure,
plot(max(dec)-dec,abs(GGev)),
hold on,
plot(dec,abs(GGmv),'--'),

figure,
semilogy(dec,1-abs(GGev).^2),
hold on,
semilogy(dec,1-abs(GGmv).^2,'--'),

figure,
semilogy(dec,abs(GGev).^2-abs(GGmv).^2), hold on
semilogy(dec,-abs(GGev).^2+abs(GGmv).^2,'.'),
semilogy(dec,1-abs(GGev(:,1)).^2),

figure,
plot(Kpf,abs(GGek)),
hold on,
plot(Kpf,abs(GGmk),'--'),

figure,
semilogy(Kpf,1-abs(GGek).^2),
hold on,
semilogy(Kpf,1-abs(GGmk).^2,'--'),
axis([0 max(Kpf) 8e-4 .8])

for k=1:8
figure,
k
semilogy(Kpf,1-abs(GGek(k,:)).^2),
hold on,
semilogy(Kpf,1-abs(GGmk(k,:)).^2,'--'),
axis([0 max(Kpf) 8e-4 .8])
pausak
end




