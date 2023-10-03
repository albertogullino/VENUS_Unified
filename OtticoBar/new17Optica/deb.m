if itetm>0
 iLP1=0;
else
 iLP1=iLP;
end


 if igamveb==1  %coeff. di riflessione sotto
   [GGe2,GGm2,TTe2,TTm2]=gaemms(KK,freq,lambda,Lvbr,nvbr,...
                   Lbb,nbb,nstratid,rfd,rr,iLP,Luvb,nuvb);
 else
  GGe2=ones(size(KK))*GGbext;
  GGm2=GGe2;
  TTe2=GGe2;
  TTm2=GGe2;
 end

 if igamveu==1  %coeff. di riflessione sopra
   [GGe1,GGm1,TTe1,TTm1]=gaemms(KK,freq,lambda,Lvtr,nvtr,...
                     Lbt,nbt,nstratiu,rfu,rr,iLP,Luv,nuv);
 else
  GGe1=ones(size(KK))*GGuext;
  GGm1=GGe1;
  TTe1=GGe1;
  TTm1=GGe1;
 end

 npk=length(KK);
 mac=1;
 mal=rr/ra;
 Zve=mac./sqrt(1-KK.^2*mac^2);
 Zvm=mac.*sqrt(1-KK.^2*mac^2);
 be0=sqrt(1-KK.^2*mac^2);




%%% B: creo dei vettori di dimensioni compatibil con il problema finale
%%%    e quindi devo replicarli nubes volte


 nubes=nubesu-nubesi;
 mbv=[-1+nubesi:nubesu+1];
 mbvz=[nubesi:nubesu];
 KKv=[];
 ZEv=[];
 ZMv=[];
 Ge1=[];
 Gm1=[];
 Ge2=[];
 Gm2=[];
 Te1=[];
 Tm1=[];
 Te2=[];
 Tm2=[];
 bev=[];
% beav=[];

 for imu=1:pasnu:nubes+1
   Ge1=[Ge1; GGe1];
   Gm1=[Gm1; GGm1];
   Ge2=[Ge2; GGe2];
   Gm2=[Gm2; GGm2];
   Te1=[Te1; TTe1];
   Tm1=[Tm1; TTm1];
   Te2=[Te2; TTe2];
   Tm2=[Tm2; TTm2];
   KKv=[KKv; KK];
   bev=[bev; be0];
%  beav=[beav; ba0];
   ZEv=[ZEv; Zve];
   ZMv=[ZMv; Zvm];
 end

 KKva=[];       %(vettore che serve a definire la sovradiagonale -1 blocco)
 nustart=2*(nube/2-fix(nube/2));
 for imu=nustart:pasnu:nubesu-2
   KKva=[KKva; KK];
 end


%%% C: Unendo TE e TM, ho i vettori totali delle dimensioni della matrice M
%%     e quindi degli autovettori a cui associo i Gamma

 Gas=[Ge1; Gm1]; %(questo verra' poi ridefinito in reluatgu)
 Gad=[Ge2; Gm2];

 Trs=[Te1; Tm1];
 Trd=[Te2; Tm2];

 be=[bev; bev];




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PARTE 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Costruzione della matrice K. Il contributo dell'anisotropia in fondo al file:
%% even Ideltaap, odd Ideltaam.

ifaccsemat=1;

if ifaccsemat==0
 if istrumix==0
  matruas
 else
  mat_mix
 end
 pes=pes(Pusa);


   if ifp>1
      disp(' ho le matrici')
      pausak
   end
end

 Gas=Gas(Pus);
 Gad=Gad(Pus);
 Trs=Trs(Pus);
 Trd=Trd(Pus);
 be=be(Pusa);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PARTE 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Su un ciclo esterno di prima even e poi odd vado in relautgu a definire
% il problema gli autovalori nei 2 casi e a risolverlo



for pol=pvet
disp(' ')
disp([' &&&&&&& Pol= ',num2str(pol)])
   if pol==1
      Iacc=Iaccp;
      KA=KAp;
      Kos=Kosp;

      if iztm==1
       KAz=KAzp;
       Kosz=Koszp;
       if iany==2
        Kosan=Kosanp;
       elseif iany==1
        Kosan=Ideltaap;
       end

       if ianys==2
        Kosan1=Kosanp;
       elseif ianys==1
        Kosan1=Ideltaap;
       end
      end

      if istrumix==1
       KTemp=KTem_p;
       if iztm==1
        KTempz=KTem_z;
       end
      else
       KTempz=0;
      end

   elseif pol==-1
      Iacc=Iaccm;
      KA=KAm;
      Kos=Kosm;

      if iztm==1
        KAz=KAzm;
        Kosz=Koszm;
        if iany==2
         Kosan=Kosanm;
        elseif iany==1
         Kosan=Ideltaam;
        end
        if ianys==2
         Kosan1=Kosanm;
        elseif ianys==1
         Kosan1=Ideltaam;
        end
      end

      if istrumix==1
       KTemp=KTem_m;
       if iztm==1
        KTempz=KTem_z;
       end
      else
       KTempz=0;
      end
   end

% relazione di dispersione

   relautan

%    disp(' relautan'), keyboard

 if pol==-1
   Gvav1=Gvav;
   alphavv1=alphavv;
   if istrumix==0
    Anvet1=Anvet;
    Ancav1=Ancav;
    if exist('Anmet')
       Anmet1=Anmet;
    end
   else
    Anz1=Anz;
    Anzf1=Anzf;
   end
 else
   Gvav2=Gvav;
   alphavv2=alphavv;
   if istrumix==0
    Anvet2=Anvet;
    Ancav2=Ancav;
    if exist('Anmet')
       Anmet2=Anmet;
    end
   else
    Anz2=Anz;
    Anzf2=Anzf;
   end

 end

%      disp(p), pausak
end  % for



clear Tstorb Tstord Tstort
