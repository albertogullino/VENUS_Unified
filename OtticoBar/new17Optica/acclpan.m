%'passo', keyboard
%iem_like=0; %0, media, 1,TE, 2,TM, -1, val in 0
%iem_like=1; %0, media, 1,TE, 2,TM, -1, val in 0
%iem_like=2; %0, media, 1,TE, 2,TM, -1, val in 0
itetm=0;

iem_like=itetm; %0, media, 1,TE, 2,TM, -1, val in 0
%' acclpan', pausak, keyboard
% [GGe1]=gaemLP(KK,pis1,pis2,pic,nstratiu,rfu);
% [GGe2]=gaemLP(KK,pis1,pis2,pic,nstratid,rfd);
if itetm>0
 iLP1=0;
else
 iLP1=1;
end
iLP1=iLP;
%keyboard
flp1=1-iLP1;
flp1=0;
iLP2=iLP1;
kv=KK;
if ilossk(1)==1
fapes=(1-1./(1+(kv/kl).^exp_los))*perdk*(-j);
else
%fapes=(1-1./(1+(kv/kl).^exp_los))*0;
fapes=zeros(size(kv));
end
if ilossk(2)==1
fapeu=(1-1./(1+(kv/kl).^exp_los))*perdk*(-j);
else
%fapeu=(1-1./(1+(kv/kl).^exp_los))*0;
fapeu=zeros(size(kv));
end


% igamveb igamveu GGbext GGuext
%'veb', keyboard
if igamveb==1
   [GGe2,GGm2,TTe2,TTm2]=gaperd(KK,freq,lambda,Lvbr,nvbr,...
                   Lbb,nbb,nstratid,rfd,rr,iLP2,Luvb,nuvb,rr,fapes);
  if iem_like==0
   GGe2=(GGe2+GGm2)/2;
   TTe2=(TTe2+TTm2)/2;
  elseif iem_like==1
% TE
   GGe2=GGe2;
   TTe2=TTe2;
% TM
  elseif iem_like==2
   GGe2=GGm2;
   TTe2=TTm2;
  elseif iem_like==-1
   GGe2=GGm2(1)*ones(size(GGe2));
   TTe2=TTm2(1)*ones(size(GGe2));
  end
%   'gaemms', keyboard
%  GGe2
else
  GGe2=ones(size(KK))*GGbext;
  GGm2=GGe2;
  TTe2=GGe2;
  TTm2=GGe2;
end
%' passp qui ', keyboard

if igamveu==1
   [GGe1,GGm1,TTe1,TTm1]=gaperd(KK,freq,lambda,Lvtr,nvtr,...
                     Lbt,nbt,nstratiu,rfu,rr,iLP2,Luv,nuv,rr,fapeu);
  if iem_like==0
   fi=find(abs(abs(GGe1)-1)<=1e-6);
   GGe1d=(GGe1+GGm1)/2;
   TTe1d=(TTe1+TTm1)/2;
   GGe1d(fi)=GGe1(fi);
   TTe1d(fi)=TTe1(fi);
   GGe1=GGe1d;
   TTe1=TTe1d;
  elseif iem_like==1
% TE
   GGe1=GGe1;
   TTe1=TTe1;
  elseif iem_like==2
% TM
   GGe1=GGm1;
  elseif iem_like==-1
   GGe1=GGm1(1)*ones(size(GGe2));
   TTe1=TTm1(1)*ones(size(GGe2));
  end

else
  GGe1=ones(size(KK))*GGuext;
  GGm1=GGe1;
  TTe1=GGe1;
  TTm1=GGe1;
end


%  [GGe2,GGm2]=gaemms(KK,freq,lambda,Lvbr,nvbr,Lbb,nbb,nstratid,rfd,rr,iLP1,...
%              Luvb,nuvb);
%  [GGe1,GGm1]=gaemms(KK,freq,lambda,Lvtr,nvtr,Lbt,nbt,nstratiu,rfu,rr,iLP1,...
%              Luv,nuv);


  if ifp==0
   disp(' verifica gamma ')
   pausak
  end
 if ifp~=-4 & ifr==1
  figure, plot(KK,real(GGe1),KK,imag(GGe1)), hold on
   plot(KK,abs(GGe1),'y.-'),
   title(' UP: Real (yellow), Imag (magenta), Abs (dots) ' ),
   xlabel(' Refl. Coefficient ')
   pausak
  figure, plot(KK,real(GGe2),KK,imag(GGe2)), hold on
   plot(KK,abs(GGe2),'y.-'),
   title(' LOW: Real (yellow), Imag (magenta), Abs (dots) ' ),
   xlabel(' Refl. Coefficient ')
   pausak
 end

%'here'
%keyboard

 npk=length(KK);
% mac=rr/real(r0);
 mac=1;
 mal=rr/ra;
 be0=conj(sqrt(1-KK.^2));
 if itetm<=0
  Zve=mac*be0./be0;
 elseif itetm==1
  Zve=mac./be0;
 elseif itetm==2
  Zve=mac.*be0;
  GGe1=GGm1;
  GGe2=GGm2;
 end

 nubes=nubesu-nubesi;
 mbv=[nubesi:nubesu];

 KKv=[];
 ZEv=[];
 Ge1=[];
 Ge2=[];
 Te1=[];
 Te2=[];
 bev=[];
 beav=[];

 for imu=1:pasnu:nubes+1
  KKv=[KKv; KK];
  ZEv=[ZEv; Zve];
  Ge1=[Ge1; GGe1];
  Ge2=[Ge2; GGe2];
  Te1=[Te1; TTe1];
  Te2=[Te2; TTe2];
  bev =[bev ;be0];
%  beav=[beav ;ba0];
 end

  if ipolar==0
    irepm=2;
  else
    irepm=1;
  end
  ZMv=ZEv;
  Gas=repmat(Ge1,irepm,1);
  Gad=repmat(Ge2,irepm,1);
  Trs=repmat(Te1,irepm,1);
  Trd=repmat(Te2,irepm,1);
  be=repmat(bev,irepm,1);

if ifaccsemat==0

 if istrumix==0
  matrlpsh
 else
%  mat_mix0
%  ' dopo matmix0', keyboard
  mat_mix
 end
end

 if irid_bas==1
%  Gas=Gas(Pus);
%  Gad=Gad(Pus);
%  Trs=Trs(Pus);
%  Trd=Trd(Pus);
%  if iredmat==0
%   Pusd=Pus;
%  else
%   Pusd=Pusas;
%  end
if iredmat==0
 Pusd=Pus;
 Pusd1=Pusas;
else
 Pusd=Pusas;
 Pusd1=Pusas;
end

  Gas=Gas(Pusd);
  Gad=Gad(Pusd);
  Trs=Trs(Pusd);
  Trd=Trd(Pusd);
%' be ', pausak

  be=be(Pusd1);
  if (freq==fre_camp(1) | iKexact==1) & icont_it==1
    pes=pes(Pusd1);
  end
 end


if ipolar==0
% Gas=[Gas; Gas];
% Gad=[Gad; Gad];
% Trs=[Trs; Trs];
% Trd=[Trd; Trd];
% be=[be; be];
 if ifr==1
%  Pust=[Pust Pust+length(Pust)];
%  pes=[pes; pes];
  Pus0=Pus;
 end
end
if ~exist('ifiez')
 ifiez=0;
end
%' acclapna', keyboard
icpo=1;
for pol=pvet

 if pol==1 | pol==0
  KA=KAp;
  Iacc=Iaccp;
  PUrD=PUriP;
  Kos=Kosp;
  KTemp=KTem_p;
  if iztm==1 & Tdis~=0
   KTempz=KTemz_p;
  else
   KTempz=0;
  end
  if itetm~=2
   KAz=0;
   Kosz=0;
  else
   KAz=KAzp;
   Kosz=Koszp;
  end
%       if ianti_gui==1
%        if N==0
%         KA_ag=KAp;
%        else
%         KA_ag=KAp_ag;
%        end
%       end
       if ianti_gui==1
%        if N==0
%         KA_ag=KAp;
%        else
         KA_ag=KAp_ag;
%        end
       end

 elseif pol==-1
%  disp(' attenzione !!!! : non ancora sistemata la parte meno (accret, matru ...)')
%  keyboard
  Iacc=Iaccm;
  PUrD=PUriM;
  KA=KAm;
%       if ianti_gui==1
%        if N==0
%         KA_ag=KAm;
%        else
%         KA_ag=KAm_ag;
%        end
%       end
       if ianti_gui==1
%        if N==0
%         KA_ag=KAm;
%        else
         KA_ag=KAm_ag;
%        end
       end
  Kos=Kosm;
  KTemp=KTem_m;
  if iztm==1 & Tdis~=0
   KTempz=KTemz_m;
  else
   KTempz=0;
  end
  if itetm~=2
   KAz=0;
   Kosz=0;
  else
   KAz=KAzm;
   Kosz=Koszm;
  end
 end

 relautan         %re_new

% if p==-1
%  Anvet1=Anvet;
%  Gvav1=Gvav;
%  alphavv1=alphavv;
%%  Ancav1=Ancav;
%  if exist('Anmet')
%   Anmet1=Anmet;
%  end
% end

 if ifiez==1
  return
 end
% if pol==-1
 if icpo==1 | pol==0
   Gvav1=Gvav;
   alphavv1=alphavv;
   if istrumix==0
    Anvet1=Anvet;
    Ancav1=Ancav;
    if exist('Anmet')
       Anmet1=Anmet;
    end
   else
%    'Anz'
%    keyboard
    Anz1=Anz;
    Anzf1=Anzf;
   end
%  elseif pol==1 | pol==0
  elseif icpo==2
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
 icpo=icpo+1;


end

clear Tstord
if ifast==0
 clear Tstorb  Tstort
end

