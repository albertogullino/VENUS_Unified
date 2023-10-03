%load qua
%load ott
%load pro
close all

KK    = Ppol.KK    ;
KKt   = Ppol.KKt   ;
numodi= Ppol.numodi;
pasnu = Ppol.pasnu ;
lbv   = Ppol.lbv   ;
kcav0 = Ppol.kcav0 ;
x     = Ppol.x     ;
fian  = Ppol.fian  ;
Nx    = Ppol.Nx    ;
mbv   = Ppol.mbv   ;
Azvet1= Ppol.Az1;
Azvet2= Ppol.Az2;
Gvet1 = Ppol.G1 ;
Gvet2 = Ppol.G2 ;
alvet1 = Ppol.A1 ;
alvet2 = Ppol.A2 ;


imod=0;
%for pola=[-1 1]
for pola=pola0

  if pola==-1
   Azvet=Azvet1;
%   Azvetf=Azvetf1;
   Gvet=Gvet1;
   alvet=alvet1;
  else
   Azvet=Azvet2;
%   Azvetf=Azvetf2;
   Gvet=Gvet2;
   alvet=alvet2;
  end
  ifps=ifp;
  if ifp>=-2
   ifp=-1;
  end
  isav=0;
  camv_new
  if ifp>=-3 | ifp==-10 | ifp==-4
   fsaf=figure,
   vg=0.5;
   puf=1:length(fso);
   subplot(211), plot(fou(puf,:)',aou(puf,:)',fso,fso*0,'wo'), grid;
   dista=(max(gou(puf,:))-min(gou(puf,:)))/min(gou(puf,:));
   title(' polar=2 ')
   if dista>10
    subplot(212), semilogy(fou(puf,:)',gou(puf,:)'/vg,fso,gso/vg,'wo');
   else
    subplot(212), plot(fou(puf,:)',gou(puf,:)'/vg,fso,gso/vg,'wo');
   end
   raf=[nomeFs(1:end-4),num2str(iLP),Ev_or_Od,'P2'];
   eval(['hgsave(',num2str(fsaf),',''',raf,''');']);
  end

end

