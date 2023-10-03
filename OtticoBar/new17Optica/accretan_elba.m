%' accretan', keyboard
iclo=1;
iem_like=10;
%icmt=0;
if itetm>0
 iLP1=0;
else
 iLP1=iLP;
end
kv=KK;
if ilossk(1)==1
fapes=(1-1./(1+(kv/kl).^exp_los))*perdk*(-j);
else
fapes=(1-1./(1+(kv/kl).^exp_los))*0;
end
if ilossk(2)==1
fapeu=(1-1./(1+(kv/kl).^exp_los))*perdk*(-j);
else
fapeu=(1-1./(1+(kv/kl).^exp_los))*0;
end

%' di qui passo ', keyboard


 if igamveb==1  %coeff. di riflessione sotto
   [GGe2,GGm2,TTe2,TTm2]=gaperd(KK,freq,lambda,Lvbr,nvbr,...
                   Lbb,nbb,nstratid,rfd,rr,iLP,Luvb,nuvb,rr,fapes);


%   lambdac=lambda*(1-freq);
%   [GGe2,GGm2,TTe2,TTm2]=gam_cmu(KK,lambdac,Lvbr,nvbr,...
%                     Lbb,nbb,nstratiu,Luvb,nuvb,rfd,rr,rr);

 else
  GGe2=ones(size(KK))*GGbext;
  GGm2=GGe2;
  TTe2=GGe2;
  TTm2=GGe2;
 end

 if igamveu==1  %coeff. di riflessione sopra

   [GGe1,GGm1,TTe1,TTm1]=gaperd(KK,freq,lambda,Lvtr,nvtr,...
                     Lbt,nbt,nstratiu,rfu,rr,iLP,Luv,nuv,rr,fapeu);

%   lambdac=lambda*(1-freq);
%   [GGe1,GGm1,TTe1,TTm1]=gam_cmu(KK,lambdac,Lvtr,nvtr,...
%                     Lbt,nbt,nstratiu,Luv,nuv,rfu,rr,rr);

 else
  GGe1=ones(size(KK))*GGuext;
  GGm1=GGe1;
  TTe1=GGe1;
  TTm1=GGe1;
 end
%'verifica', keyboard
%shafi=shavet(2:end-1,1);
%fili=find(shafi~=0 & aitot==0);

 if ifp~=-4 & ifr==1
  figure, plot(KK,real(GGe1),KK,imag(GGe1)), hold on
   plot(KK,real(GGm1),'--',KK,imag(GGm1),'--'),
   plot(KK,abs(GGm1),'m.-',KK,abs(GGe1),'y.-'),
   title(' TE and TM(--) UP: Real (yellow) and Imag (magenta) ' ),
   xlabel(' Refl. Coefficient ')
   pausak
  figure, plot(KK,real(GGe2),KK,imag(GGe2)), hold on
   plot(KK,real(GGm2),'--',KK,imag(GGm2),'--'),
   plot(KK,abs(GGm2),'m.-',KK,abs(GGe2),'y.-'),
   title(' TE and TM(--) LOW: Real (yellow) and Imag (magenta) ' ),
   xlabel(' Refl. Coefficient ')
   pausak
   if iclo==1
   close(1:2)
   end
 end
% 'gaemms', keyboard

 npk=length(KK);
 mac=1;
 mal=rr/ra;
 be0=conj(sqrt(1-KK.^2*mac^2));
 Zve=mac./be0;
 Zvm=mac.*be0;




%%% B: creo dei vettori di dimensioni compatibil con il problema finale
%%%    e quindi devo replicarli nubes volte


 nubes=nubesu-nubesi;
 mbv=[-1+nubesi:nubesu+1];
% 'mbv', keyboard
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
% ZEv=real(ZEv);
% ZMv=real(ZMv);


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

 if ipolar==0
    irepm=2;
  else
    irepm=1;
  end
  Gas=repmat(Gas,irepm,1);
  Gad=repmat(Gad,irepm,1);
  Trs=repmat(Trs,irepm,1);
  Trd=repmat(Trd,irepm,1);
  be=repmat(be,irepm,1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PARTE 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Costruzione della matrice K. Il contributo dell'anisotropia in fondo al file:
%% even Ideltaap, odd Ideltaam.

%  '  matmix', keyboard
if ifaccsemat==0
 if istrumix==0
  matruas
 else
%  mat_mix0
%global imatm
%'mat', keyboard
%  if length(imatm)==0
   mat_mix
%   imatm=1;
%   save matmix
%  else
%   save par
%   load matmix
%   load par
%  end 
 end
 pes=pes(Pusas);


   if ifp>1
      disp(' ho le matrici')
      pausak
   end
end
if iredmat==0
 Pusd=Pus;
 Pusd1=Pus;
else
 Pusd=Pusas;
 Pusd1=Pusas;
end
%' passo', keyboard

 Gas=Gas(Pusd);
 Gad=Gad(Pusd);
 Trs=Trs(Pusd);
 Trd=Trd(Pusd);

 be=be(Pusd1);

%if ifp==-10 & ifr==1
% ' salvo tutto '
% save par
% keyboard
%end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PARTE 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Su un ciclo esterno di prima even e poi odd vado in relautgu a definire
% il problema gli autovalori nei 2 casi e a risolverlo


icpo=1;

if ~exist('pvet00')
 pvet00=pvet;
end
pcounter=0;
for pol=pvet00
pcounter=pcounter+1;

if ifp~=-4
disp(' ')
disp([' &&&&&&& Pol= ',num2str(pol)])
end
   if pol==1 | pol==0
      Iacc=Iaccp;
      PUrD=PUriP;
      KA=KAp;
      KAN=KANp;
      Kos=Kosp;
%      ' qui ver: pol=1 ', keyboard
       if ianti_gui==1
        if N==0
         KA_ag=KAp;
         KAz_ag=KAzp;
        else
         KA_ag=KAp_ag;
         KAz_ag=KAzp_ag;
        end
       end

%      ' qui ver dopo', keyboard
      if iztm==1
       KAz=KAzp;
       Kosz=Koszp;
       KANz=KANzp;
%       if ianti_gui==1
%        KAz=KAzp;
%        Kosz=Koszp;
%       end
      end

       if iany==2
        Kosan=Kosanp;
       elseif iany==1
        Kosan=Ideltaap;
       end
       if exist('i_grap')
        if i_grap==2
         Kan_lim=Kanpp;
        elseif i_grap==1
         Kanr=Ideltaap;
        end
       end

       if ianys==2
        Kosan1=Kosanp;
       elseif ianys==1
        Kosan1=Ideltaap;
       end

      if istrumix==1
       KTemp=KTem_p;
       if iztm==1
        KTempz=KTep_z;
       end
      else
       KTempz=0;
      end

   elseif pol==-1
      Iacc=Iaccm;
      KAN=KANm;
      PUrD=PUriM;
      KA=KAm;
      Kos=Kosm;
%      ' qui ver: pol=-1 ', keyboard
       if ianti_gui==1
        if N==0
         KA_ag=KAm;
         KAz_ag=KAzm;
        else
         KA_ag=KAm_ag;
         KAz_ag=KAzm_ag;
        end
       end


      if iztm==1
       KAz=KAzm;
       Kosz=Koszm;
       KANz=KANzm;
%       if ianti_gui==1
%        KAz=KAzm;
%        Kosz=Koszm;
%       end
      end

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
       if exist('i_grap')
        if i_grap==2
         Kan_lim=Kanpm;
        elseif i_grap==1
         Kanr=Ideltaam;
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

%    disp(' relautan in accretan'), pausak
%    pausak
   relautan


% if pol==-1
 if icpo==1
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
%   'set Anz1', pausak
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
%   'set Anz2', pausak

 end
 icpo=icpo+1;

%      disp(p), pausak
end  % for



clear Tstord
if ifast==0
 clear Tstorb  Tstort
end
