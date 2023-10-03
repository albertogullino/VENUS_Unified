%ad
iade=0;    % fa come in programma
clear Im
ifp=-10;
icam_fr=input(' campi per ogni frequenza ? [0/1] ');
if length(icam_fr)==0
 icam_fr=0;
end
ionly1=input(' solo una soluzione ? [0/1] ');
if length(ionly1)==0
 ionly1=0;
else
 imod=input(' modo 1 o 2 ? [1/2] ');
 if imod==1
  ipolar=1;
 else
  ipolar=-1;
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

  if ipolar==2
   pola=pvet(1);
   Azvet=Azvet1;
   Azvetf=Azvetf1;
   Gvet=Gvet1;
   alvet=alvet1;
%   ' qui ipolar', keyboard

  else
   pola=ipolar;
   if ipolar==1 | ipolar==0
    Azvet=Azvet1;
    Azvetf=Azvetf1;
    Gvet=Gvet1;
    alvet=alvet1;
   else
%   elseif ipolar==1 | ipolar==0
    Azvet=Azvet2;
    Azvetf=Azvetf2;
    Gvet=Gvet2;
    alvet=alvet2;
   end
  end

  imod=0;
  isalva=1;
  ' rima ', keyboard
 if iade==0
  camv_new
 else
  ca_de
 end
%  'ferma'
%  keyboard
   lco=lambda*1000;
saou=1
  if (ifp>=-3 | ifp==-10 | ifp==-4) & saou>0
   vg=0.5;
   fsaf=figure,
   puf=1:length(fso);

   zeromo=ones(size(fso))*alpha_th;
%c   subplot(211), plot(lco*fou(:,puf),aou(:,puf),lco*fso,zeromo,'wo'), grid;
   subplot(211), plot(lco*fou(puf,:)',aou(puf,:)',lco*fso,zeromo,'wo'), grid;
   title(' polar=1 ')
   dista=(max(gou(puf,:))-min(gou(puf,:)))/min(gou(puf,:));
   if dista>10
    subplot(212), semilogy(lco*fou(puf,:)',gou(puf,:)'/vg,lco*fso,gso/vg,'wo');
   else
    subplot(212), plot(lco*fou(puf,:)',gou(puf,:)'/vg,lco*fso,gso/vg,'wo');
   end

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
  if (ifp>=-3 | ifp==-10 | ifp==-4) & saou>0
   fsaf=figure,
   vg=0.5;
   puf=1:length(fso);
   zeromo=ones(size(fso))*alpha_th;
   subplot(211), plot(lco*fou(puf,:)',aou(puf,:)',lco*fso,zeromo,'wo'), grid;
   dista=(max(gou(puf,:))-min(gou(puf,:)))/min(gou(puf,:));
   title(' polar=2 ')
   if dista>10
    subplot(212), semilogy(lco*fou(puf,:)',gou(puf,:)'/vg,lco*fso,gso/vg,'wo');
   else
    subplot(212), plot(lco*fou(puf,:)',gou(puf,:)'/vg,lco*fso,gso/vg,'wo');
   end
  end

 end
 ifp=ifps;

%pol_rpo
