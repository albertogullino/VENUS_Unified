colordef black
%icampi=0
%ad
%close all
isk=1
if isk==0
if exist('R0v')
    R=R0v;
  O=ones(size(R0v));
Gm=mean(mean(G));
fi=find(G<Gm);
L_sup1=max(max(L));
L_inf1=min(min(L));  
Lm=mean(L(fi));
  figure
  subplot(211), plot(R,L*1000,'.-',R,O*L_inf*1000,'r',R,O*L_sup*1000,'g')
  hold on, plot(R,O*L_inf1*1000,'r--',R,O*L_sup1*1000,'g--',R,O*Lm*1000,'c-o')
  xlabel([' radial coord. (micron)'])
  subplot(212), plot(R,G,R,O*Gm,'r')
  pausak
end
end
clear Gsov Fsov
isavetutto=0;
%ipolar=2;
%nmasce=1;
%'postg', keyboard
ipostg=1;        %=0 salva
ipost=0;        %=0 ricalcola base
isav_Az=0;        %=0 calcolo campi z
Ps.isav_Az=isav_Az;        %=0 calcolo campi z
%ipost=0;        %=0 ricalcola base
iade=0;    % fa come in programma
%clear Im
ifp=-10;
%ifp=-4;
%icam_fr=input(' campi per ogni frequenza ? [0/1] ');
if ~exist('icam_fr')
 icam_fr=0;
end
%ionly1=input(' solo una soluzione ? [0/1] ');
ionly1=1;
if length(ionly1)==0
 ionly1=0;
else
ipolarv=ipolar;
% imod=input(' modo 1 o 2 ? [1/2] ');
 imod=1;
 ICON=input(' ICON  ? [0/1] ');
if length(ICON)==0
 ICON=0;
end
% ICON=0;
 if imod==1
  pvet=1;
  ipolar=pvet;
 else
  pvet=-1;
  ipolar=pvet;
 end
end


%clear global
%load
%ipolar=-1;
%  ifp=-10;
%  if exist('Pus0')
%   Pus=Pus0;
%  end

  ifps=ifp; if ifp>=-2; ifp=1; end
     pola=pvet(1);
     Azvet=Azvet1;
     Azvetf=Azvetf1;
     Gvet=Gvet1;
   alvet=alvet1;
if ipolarv==2 & ipolar==1
   pola=pvet(1);
   Azvet=Azvet1;
   Azvetf=Azvetf1;
   Gvet=Gvet1;
   alvet=alvet1;
elseif ipolarv==2 & ipolar==-1
   pola=pvet(1);
   Azvet=Azvet2;
   Azvetf=Azvetf2;
   Gvet=Gvet2;
   alvet=alvet2;
end

  imod=0;
  isalva=1;
  ' rima ', keyboard
 if iade==0
  camv_new1
 else
  ca_de
 end
%  'ferma'
%  keyboard
   lco=lambda*1000;
saou=1
  if (ifp>=-3 | ifp==-10 | ifp==-4)
   vg=0.5;
   vg=vgconv;
   vg1=vgconv/2;
   fsaf=figure;
   set(fsaf,'pos',[433    23   672   504])
   puf=1:size(aou,1);
   
   zeromo=ones(size(Fsov))*alpha_th;
%c   subplot(211), plot(lco*fou(:,puf),aou(:,puf),lco*fso,zeromo,'wo'), grid;
   subplot(211), plot(lco*fou(puf,:)',aou(puf,:)',lco*Fsov,zeromo,'wo'), grid;
   title(' polar=1 ')
   dista=(max(gou(puf,:))-min(gou(puf,:)))/min(gou(puf,:));
   dista=20;
   fim=find(M2v<M2_max);
      Gso=Gsosa/vg;
   fim=find(Gso<1.5*min(Gso));
   if dista>10
   Fso=lco*Fsov;
    subplot(212), semilogy(lco*fou(puf,:)',gou(puf,:)'/vg1,Fso,Gso,'wo',Fso(fim),Gso(fim),'cs');
   else
    subplot(212), plot(lco*fou(puf,:)',gou(puf,:)'/vg1,lco*Fsov,Gsosa/vg,'wo');
   end

   grid

  end
if ifp==-10  
for k=fim, figure(k+1), pausak, end
end
pausak
%keyboard
% icasa=1;
fig=find(gsov>0);
figure, plot(M2v(fig),gsov(fig)/vgconv,'ro'), pausak
dlam=lambda./(1+fsov(fig));
dlam=lambda*fsov*1000;
Gcm=2*gsov/vgconv;
%figure, plot(M2v(fig),dlam,'ro'), pausak
figure, semilogy(dlam(fig),gsov(fig)/vgconv,'go'), grid, pausak
figure, plot(dlam(fig),M2v(fig),'co'), grid, pausak

[du,iso]=sort(dlam);
[2*gsov(iso)/vgconv M2v(iso)' dlam(iso) [1:length(iso)]']

 if ipolar==2

   ' qui ipolar secondo', keyboard
  Azvet=Azvet2;
  Azvetf=Azvetf2;
  Gvet=Gvet2;
  alvet=alvet2;
  pola=pvet(2);
  isalva=1;
%  camv_pro
%  camv_new
 if iade==0
  camv_new1
 else
  ca_de
 end

%  Ppol.Ex2=Pol.Ex;
%  Ppol.Ey2=Pol.Ey;

%  if ifp>=-3 | ifp==-10 | ifp==-4
  if (ifp>=-3 | ifp==-10 | ifp==-4) 
   figure(fsaf),
   vg=0.5;
   puf=1:length(fso);
   zeromo=ones(size(fso))*alpha_th;
   subplot(211),  hold on
   plot(lco*fou(puf,:)',aou(puf,:)','--',lco*fso,zeromo,'ws'), grid;
   dista=(max(gou(puf,:))-min(gou(puf,:)))/min(gou(puf,:));
   dista=20;
   if dista>10
%    subplot(212), semilogy(lco*fou(puf,:)',gou(puf,:)'/vg,lco*fso,gso/vg,'wo');
%   else
    subplot(212), hold on, plot(lco*fou(puf,:)',gou(puf,:)'/vg,'--',lco*fso,gso/vg,'ws');
   end
  end

 end
 ifp=ifps;

%pol_rpo
