ifpsave=ifp;
ifp=-4;

if iraff==1
    pis1=pis10*(1+freq);
    pis2=pis20*(1+freq);
    pic=pic0*(1+freq);
    kcav=kcav0*(1+freq);
    lax=(lambda+Dla)/(1+fris);
    acc0_i
%    eval(['load ',nfsave])

    fim=find(Grlam>0 & Grlam<10*gg0);
    Corre=abs(vet(:,fim)'*ADoutis);

    if ifpsave>=-1
     disp(' G alfa Corre ')
     [Grlam(fim) allam(fim) Corre]
    end
    [du,fsold1]=min(abs(allam(fim)));

    [du,fsold]=max(Corre);
    if fsold~=fsold1
     disp(' diversi')
%     pausak
    end
    fsol=fim(fsold);

    Ar=vet(:,fsol);   % autov all'interfaccia
    Anout=Ar;
    glosout=2*Grlam(fsol);
    [2*gg0 glosout]

      Apr=[Gas.*Ar; Ar];
      CCd=Apr;

      ics=1;
      Ab=CCd(pu2);
      Af=CCd(pu1);
      Afz(:,ics)=Ab.*Trs;
      Afzf(:,ics)=Af.*Trs;

      CCd=Tt*Apr;
      ics=ics+1;
      Ab=CCd(pu2);
      Af=CCd(pu1);
      Afz(:,ics)=Ab;
      Afzf(:,ics)=Af;

      Lam=Lamv(iN);
      Ttotal=(Mn+Lam*Mi);
      CCd=Ttotal*Apr;
      ics=ics+1;
      Ab=CCd(pu2);
      Af=CCd(pu1);
      Afz(:,ics)=Ab.*Trd;
      Afzf(:,ics)=Af.*Trd;

      sdu=size(Tmef);
      sTm=sdu(1:2);

%  sopra
      if icsfit>0
       for kf=1:icsfit
        ics=ics+1;
        CCd=reshape(Tmef(:,:,kf),sTm)*Apr;
        Ab=CCd(pu2);
        Af=CCd(pu1);
        Afz(:,ics)=Ab;
        Afzf(:,ics)=Af;
       end
      end

%  sotto
      if icsfib>0
       TQW=(Mon+Lam*Moi);
       AprQW=TQW*Tt*Apr;
       for kf=1:icsfib
        ics=ics+1;
        CCd=reshape(Tmefb(:,:,kf),sTm)*AprQW;
        Ab=CCd(pu2);
        Af=CCd(pu1);
        Afz(:,ics)=Ab;
        Afzf(:,ics)=Af;
       end
      end

%disp(' fine autov '), keyboard
 clear Tb Tt Moi Mon Tmefb Tmef Mn Mi
else


   Afz=reshape(Azvet(:,pou(nso,fiAz),:,fiAz),sz([1 3]));
   Afzf=reshape(Azvetf(:,pou(nso,fiAz),:,fiAz),sz([1 3]));
   Anout=Afz(:,1)+Afzf(:,1);
   AnQW=Afz(:,2)+Afzf(:,2);
   ics=sz(3);


end

 if icampi==2
  icsmax=ics;
 else
  icsmax=2;
 end


 for kmo=1:icsmax
   Andu=Afz(:,kmo)+Afzf(:,kmo);

%  plot FF

   if kmo==1 & iFF==1
%     Anduf=Afzf(:,kmo);
%     AnFF=Andu+Anduf;
     AnFF=Andu;
     plo_ff
   end

   cam_val;
   if iskim==1, break, end
   if kmo==1
    E2x=E2xd;
    E2y=E2yd;
ifp=ifpsave;
    return
    mod_cha
    disp(' salva'), keyboard
   end
   if iskim==1, break, end
   if i2D==3
    Etb(:,:,kmo)=Etot;
    Etx(:,:,kmo)=E2xd;
    Ety(:,:,kmo)=E2yd;
   else
    Etb(:,kmo)=Etot;
    Etx(:,kmo)=E2xd;
    Ety(:,kmo)=E2yd;
   end
   if kmo==2
      Andu=Afz(:,kmo)+Afzf(:,kmo);
      cam_val;
      Etcsum=Etot;
   end

   if idyn==1
    Andu=Afzf(:,kmo);
    cam_val;
    if i2D==3
     Etf(:,:,kmo)=Etot;
    else
     Etf(:,kmo)=Etot;
    end
   end
 end


%eval(['!del ',nfsave,'.mat'])
