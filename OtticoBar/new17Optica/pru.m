clear 
load fui
ifp=-10
gg0=GStim(nso);
ze=FStim(nso);

   if gg0/vg<GMA
    nsc=nsc+1;
    ns=[ ns nso];
    if iftot(nso)==2
     ipufa=1000;
     fiAz=length(Fint);
    else
     ipufa=-1000;
     fiAz=1;
    end 
    ipu=[ipu ipufa];
    fso=[fso ze];
    gso=[gso gg0];
    aso=[aso 0];
    tso=[tso -1];
   end
   sA=length(puA);
   An0=reshape(Anu(:,pou(nso,fiAz),fiAz),sA(1),1);
   Anso=[Anso An0];
   
%' capeoi', keyboard   

   if icampi>=1

%    fieval
%   disp('camdu in diss_nst 1'), keyboard
%ifp=-10
    lep=size(pou);
    Fint=Fint(1:lep(2));
    iLP=iLP1;
    fie_new
    iLP=iLPr;

%    'ferma', keyboard
    if iLP==1
     rtetm=0.5;
     nuazi=0;
     polca=0;
     polratio=0;
     mrad=0;
    end
    M2l=[M2l M2];
    rtetmv=[rtetmv rtetm];
    maziv=[maziv nuazi];
    mradv=[mradv mrad];
    polcav=[polcav polca];
    polrat=[polrat polratio];
%        ' cont diss', keyboard
 %  disp('polcav'), pausak
   end
