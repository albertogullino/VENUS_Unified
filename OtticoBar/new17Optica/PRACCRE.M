%' accretan',
icmt=1;
if itetm>0
 iLP1=0;
else
 iLP1=iLP;
end

%' di qui passo ', keyboard


 if igamveb==1  %coeff. di riflessione sotto
  if icmt==0
   [GGe2,GGm2,TTe2,TTm2]=gaemms(KK,freq,lambda,Lvbr,nvbr,...
                   Lbb,nbb,nstratid,rfd,rr,iLP,Luvb,nuvb);
  else
   lambdac=lambda*(1-freq);
%   [GGe2,GGm2,TTe2,TTm2,Tpe2,Tpm2]=gam_cmu(KK,lambdac,Lvbr,nvbr,...
   [GGe2,GGm2,TTe2,TTm2]=gam_cmu(KK,lambdac,Lvbr,nvbr,...
                     Lbb,nbb,nstratiu,Luvb,nuvb,rfd,rr,rr);
  end
%   [GGe2,GGm2,TTe2,TTm2]=gaemms(KK/3,freq,lambda,Lvbr,nvbr,...
%                   Lbb,nbb,nstratid,rfd,rr,iLP,Luvb,nuvb);

 else
  GGe2=ones(size(KK))*GGbext;
  GGm2=GGe2;
  TTe2=GGe2;
  TTm2=GGe2;
 end

 if igamveu==1  %coeff. di riflessione sopra
  if icmt==0
   [GGe1,GGm1,TTe1,TTm1]=gaemms(KK,freq,lambda,Lvtr,nvtr,...
                     Lbt,nbt,nstratiu,rfu,rr,iLP,Luv,nuv);
  else
   lambdac=lambda*(1-freq);
%   [GGe1,GGm1,TTe1,TTm1,Tpe1,Tpm1]=gam_cmu(KK,lambdac,Lvtr,nvtr,...
   [GGe1,GGm1,TTe1,TTm1]=gam_cmu(KK,lambdac,Lvtr,nvtr,...
                     Lbt,nbt,nstratiu,Luv,nuv,nuv(end),rr,rr);
  end
%   ' prova ref '
%   keyboard
%   pga

 else
  GGe1=ones(size(KK))*GGuext;
  GGm1=GGe1;
  TTe1=GGe1;
  TTm1=GGe1;
 end
%'verifica', keyboard

 if ifp~=-4 & ifr==1
  figure, plot(KK,real(GGe1),KK,imag(GGe1)), hold on
   plot(KK,real(GGm1),'--',KK,imag(GGm1),'--'),
   title(' TE and TM(--) refl coeff upper interface' ), pausak
  figure, plot(KK,real(GGe2),KK,imag(GGe2)), hold on
   plot(KK,real(GGm2),'--',KK,imag(GGm2),'--'),
   title(' TE and TM(--) refl coeff lower interface' )
   pausak
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
