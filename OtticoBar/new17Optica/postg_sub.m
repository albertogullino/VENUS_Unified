%ad
%close all
clear Gsov Fsov
ipolar=2;
nmasce=5;
ipost=1;        %=0 ricalcola base
isav_Az=0;        %=0 calcolo campi z
Ps.isav_Az=0;        %=0 calcolo campi z
%ipost=0;        %=0 ricalcola base
iade=0;    % fa come in programma
%clear Im
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
if ifp==-10
 ICON=input(' ICON  ? [0/1] ');
else 
 ICON=0;
end 
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
%  ' rima ', keyboard
 if iade==0
  camv_new
 else
  ca_de
 end
%  'ferma'
%  keyboard
   lco=lambda*1000;
saou=1
  if (ifp>=-3 | ifp==-10 )
   vg=0.5;
   vg=vgconv;
   vg1=vgconv/2;
   fsaf=figure,
   puf=1:size(aou,1);

   zeromo=ones(size(Fsov))*alpha_th;
%c   subplot(211), plot(lco*fou(:,puf),aou(:,puf),lco*fso,zeromo,'wo'), grid;
   subplot(211), plot(lco*fou(puf,:)',aou(puf,:)',lco*Fsov,zeromo,'wo'), grid;
   title(' polar=1 ')
   dista=(max(gou(puf,:))-min(gou(puf,:)))/min(gou(puf,:));
   dista=20;
   if dista>10
    subplot(212), semilogy(lco*fou(puf,:)',gou(puf,:)'/vg1,lco*Fsov,Gsosa/vg,'wo');
   else
    subplot(212), plot(lco*fou(puf,:)',gou(puf,:)'/vg1,lco*Fsov,Gsosa/vg,'wo');
   end
   grid

  end

%keyboard
% icasa=1;

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
  camv_new
 else
  ca_de
 end

%  Ppol.Ex2=Pol.Ex;
%  Ppol.Ey2=Pol.Ey;

%  if ifp>=-3 | ifp==-10 | ifp==-4
  if (ifp>=-3 | ifp==-10 ) 
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
