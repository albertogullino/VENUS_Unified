%clear
%clear global
%load
ipolar=1;
%  ifp=-10;
%  if exist('Pus0')
%   Pus=Pus0;
%  end
global if_only_out
if_only_out=1;
ifp=-10

  ifps=ifp; if ifp>=-2; ifp=1; end

  if ipolar==2
   pola=pvet(1);
   Azvet=Azvet1;
   Azvetf=Azvetf1;
   Gvet=Gvet1;
   alvet=alvet1;

  else
   pola=ipolar;
   if ipolar==-1
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
  camv_new
%  'ferma'
%  keyboard

  if ifp>=-3 | ifp==-10 | ifp==-4
   vg=0.5;
   fsaf=figure,
   puf=1:length(fso);

   zeromo=ones(size(fso))*alpha_th;
   subplot(211), plot(fou(puf,:)',aou(puf,:)',fso,zeromo,'wo'), grid;
   title(' polar=1 ')
   dista=(max(gou(puf,:))-min(gou(puf,:)))/min(gou(puf,:));
   if dista>10
    subplot(212), semilogy(fou(puf,:)',gou(puf,:)'/vg,fso,gso/vg,'wo');
   else
    subplot(212), plot(fou(puf,:)',gou(puf,:)'/vg,fso,gso/vg,'wo');
   end

  end


% icasa=1;

 if ipolar==2

  Azvet=Azvet2;
  Azvetf=Azvetf2;
  Gvet=Gvet2;
  alvet=alvet2;
  pola=pvet(2);
  isalva=1;
  camv_new

%  Ppol.Ex2=Pol.Ex;
%  Ppol.Ey2=Pol.Ey;

  if ifp>=-3 | ifp==-10 | ifp==-4
   fsaf=figure,
   vg=0.5;
   puf=1:length(fso);
   zeromo=ones(size(fso))*alpha_th;
   subplot(211), plot(fou(puf,:)',aou(puf,:)',fso,zeromo,'wo'), grid;
   dista=(max(gou(puf,:))-min(gou(puf,:)))/min(gou(puf,:));
   title(' polar=2 ')
   if dista>10
    subplot(212), semilogy(fou(puf,:)',gou(puf,:)'/vg,fso,gso/vg,'wo');
   else
    subplot(212), plot(fou(puf,:)',gou(puf,:)'/vg,fso,gso/vg,'wo');
   end
  end

 end
 ifp=ifps;

%pol_rpo
