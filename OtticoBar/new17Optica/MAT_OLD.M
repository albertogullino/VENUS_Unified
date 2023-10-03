anyret=zeros(size(dv));
anyr.t=anyret(put,:);
anyr.b=anyret(pub,:);
global ikr
%'qui mat_mix ingresso', keyboard
ireturn=0;

%global inuo_bas iredmat nucomp
global inuo_bas nucomp
clear Kosp Kosm Koszp Koszm
%' entro mix', keyboard


%if length(inuo_bas)==0
% inuo_bas=0;
%end
inuo_bas=1;
ianmat=1;


if inuo_bas==1
 iredmat=1;
else
 iredmat=0;
 nucomp=0;
end
if iLP==1
 Pusas=1:(numodiacc+1)*length(KK);
else
 Pusas=1:(numodiacc+1)*length(KK)*2;
end

if ifr==1
 Tfas=0;
 Ffas=0;
 if ifp~=-4
 global ptim
  if length(ptim)==0
   ptim=figure;
  end
 end
end
 [Tfas,Ffas]=eltime(Tfas,Ffas,-10);


icosha=0;
clear Cug


%fis0=find(radii.a(:,1)==0);
fis0=find(radii.a==0);
sizsh=size(shavet);
dvc=repmat(dv,1,sizsh(2));
shavet(fis0)=0;

%fish=find(shavet(:,1)~=0 & dvc~=0);
%[shai,ishd]=sort(shavet(fish,1)));
fish=find(shavet~=0 & dvc~=0);
[shai,ishd]=sort(shavet(fish));
ish=fish(ishd);
%' mix ', keyboard

fis0du=find(radii.a(:,1)==0);
shavetdu=shavet(:,1);
fishdu=find(shavetdu~=0 & dv~=0);
shavetdu(fis0du)=0;
[shaidu,ishddu]=sort(shavetdu(fishdu));
ishdu=fishdu(ishddu);
ishdu=repmat(ishdu,sizsh(2),1);
%'ish', keyboard


ax=radii.b(ish);
ay=radii.a(ish);
pad=radii.c(ish);
sha_ish=shavet(ish);
%firet=find(shavet(ish)==6 & dv();
firet=find(shavet(ish)==6);

if length(firet)>=1
 for kfir=firet'
  icel=ish(kfir);
  DU=radii.array{icel};
  pad(kfir)=DU{1};

 end

%end  % prima era qui


%'quit'
%keyboard
 ary=radii.array(ish,:);
 ico=0;

 for duish=1:length(ary)
  ico=ico+1;
  du=ary{ico}{1};
  if length(du)==0
   du=0;
  end
  ty(ico,1)=du;
 end

end
%end  % prima era qui


 ity=0;
if length(find(shavet>4 & shavet~=7 & shavet~=6))>0
%if length(find(abs(shavet)>4 & shavet~=7 & shavet~=6))>0
 ishve=find(shavet>4 & shavet~=7);
 ishve=find(shavet>0 & shavet~=7);
% ishve=find(abs(shavet)>0 & shavet~=7);
% ary=radii.array(ishdu);
 ary=radii.array(ishve);
 ico=0;
% for duish=ishdu'
 for duish=ishve'
  ico=ico+1;
  ty(ico,1)=ary{ico}{1};
 end
 ity=1;
end


sait=size(aitot);
if sait(2)>3
 ity=0;
end

sis=size(shavet);
ifil=ifield;
% prima c'era
for k=2:sis(2)
 ifil=[ifil ifield];
% ish=[ish ish];
end

%fir=find(shavet==6 & dv~=0);
fir=find(shavet==6);
if ikr==0
 fir=[];
end
%' fir=find(shavet==6); ', keyboard
%' mat_mix ', pausak


if length(nucomp)==0
 nucomp=0;
end
if nucomp==0
 inume=0;
else
 inume=1;
end
if ifp==-10
' prima riretu', keyboard
end
if inume==1
 if ifr==1 & iredmat==1
% if length(fir)>0
  ri_retu
  si2sav=si2;
% end
 end

end
 if ifp~=-4 & exist('Pus0')
   Pud=Pus0;
  figure,
  subplot(211)
  plot(KKv(Pud),real(Gas(Pud)),'.',KKv(Pud),real(Gad(Pud)),'.'),
  subplot(212)
  plot(KKv(Pud),imag(Gas(Pud)),'.',KKv(Pud),imag(Gad(Pud)),'.'),
  pausak
 end

ifipu=ifil(ish);
ifipud=ifipu;

%fish=find(shavet~=0 & dvc~=0);
[shai,ishd]=sort(shavet(fish));
fisha=find(diff([0; shai])~=0);
shailoop=shai(fisha)';
shailoop=abs(shailoop);
shai=abs(shai);
for shad=shailoop
 if shad<0
  shad=-shad;
 end
 sha=abs(shad);
 if ifp~=-4
 disp([' mat_mix: sha = ',num2str(sha)])
 end
% ' saha ', keyboard
 fis=find(shai==sha);
   if sha==7
    sha=1;
   end
 if length(fis)>0
 if sha~=8
  if sha<5
   [ayi,iso]=sort(ay(fis));
  else
   [ayi,iso]=sort(ty(fis));
  end
  if length(fis)>1
   fish=find(diff([0; ayi])~=0);
   isor=fis(iso(fish));
   fim=find(ifipu(fis)<0);
%   fim=find(ifipu<0);
   if length(fim)>0
    if length(fim)==length(isor)
%     ' cosa fa'
%     keyboard
%     ifipud(isor)=ifipu(fim);
    else
     fim1=find(ifipu(isor)<0);
     if length(fim1)<length(fim)
%      disp(' problema in mat_mix')
%      keyboard
     end
    end
   end
  else
   isor=fis;
  end
%disp(' mat_mix interno')
%keyboard
  axi=ax(isor);
  ayi=ay(isor);
  pdi=pad(isor);

  if ity==1
   tyi=ty(isor);
  end
  ifii=ifipud(isor);
  axtot(1:length(axi),shad)=axi;
  aytot(1:length(axi),shad)=ayi;
  pdtot(1:length(axi),shad)=pdi;
  if ity==1
   tytot(1:length(axi),1)=tyi;
  end
%disp(' mat_mix interno'), keyboard

%   xi=axi/a0ref;
%   yi=ayi/a0ref;
%   adis=yi*aSAV;
%   bdis=xi*aSAV;

   adis=ayi*kcav;
   bdis=axi*kcav;
   avero=ayi;
   bvero=axi;
   lxi=length(adis);
  end  %sha~=8
%  ' qui sha3 ', keyboard
   if sha>3
    if sha==4
     shape='sha_rhom';
    elseif sha==5
%     shape='sha_ar';
%     shape='sha_arm';

     shape='s_array';
%     disp(' mat_mix'), shape
%     keyboard

     icoso=0;
     for kisor=isor'
     icoso=icoso+1;
      tyceiv=ary{kisor}{9};
      if tyceiv==1
       cce=ary{kisor}{3};
       sh_type=ary{kisor}{4};
       Rvetx=ary{kisor}{6};
       Rvety=ary{kisor}{5};
       Delta=ary{kisor}{7};
      else
       ncex=ary{kisor}{2};
       ncey=ary{kisor}{3};
       dxce=ary{kisor}{8};
       dyce=dxce;
       Rx=ary{kisor}{6};
       Ry=ary{kisor}{5};
       sh_ty=ary{kisor}{4};
       Del=ary{kisor}{7};
       if length(ary{kisor})==10
        Tyar=ary{kisor}{10};
       else
        Tyar=0;
       end

       if length(ary{kisor})==11
        Sh_ar=ary{kisor}{11};
       else
        Sh_ar='square';
       end
%       ' Tyar '
%       keyboard

     switch Sh_ar

      case {'square'}

       ncell=ncex*ncey;
       dcex=2*Rx+dxce;
       dcey=2*Ry+dyce;
       Rvetx=repmat(Rx,1,ncell);
       Rvety=repmat(Ry,1,ncell);
       Delta=repmat(Del,1,ncell);
       sh_type=repmat(sh_ty,1,ncell);


       fcex=([1:ncex]-ncex/2-.5)*dcex;
       fcey=([1:ncey]-ncey/2-.5)*dcey;
       Fcx=ones(size(fcey'))*fcex;
       Fcy=fcey'*ones(size(fcex));
       cced=Fcx+j*Fcy;

       cce=reshape(cced,1,ncell);

       if Tyar~=0
        ic=Tyar;
        icv=0;
        for ncc=1:ncex
         for ncr=1:ncey
          ic=ic+1;
          if is_even(ic)
           icv=icv+1;
           cceu(icv)=cced(ncc,ncr);
          end
         end
        end
        cce=cceu;
       end

      case {'pc'}

       mx=ncex;
       my=ncey;
       ncex=(mx-1)*2+1;
       ncey=my;
       if is_even((mx-1)/2)==0
        Tyar=-1;
       else
        Tyar=1;
       end
       if Tyar==-1
        Tyar=2;
       end

       ncell=ncex*ncey;
       dcex=dxce/2;
       dcey=dyce*sqrt(3)/2;

       Rad_max=dcex*(mx-1);

       Rvetx=repmat(Rx,1,ncell);
       Rvety=repmat(Ry,1,ncell);
       Delta=repmat(Del,1,ncell);
       sh_type=repmat(sh_ty,1,ncell);

       fcex=([1:ncex]-ncex/2-.5)*dcex;
       fcey=([1:ncey]-ncey/2-.5)*dcey;
       Fcx=ones(size(fcey'))*fcex;
       Fcy=fcey'*ones(size(fcex));
       cced=Fcx+j*Fcy;

       if Tyar~=0
        ic=Tyar;
        icv=0;
        for ncc=1:ncex
         for ncr=1:ncey
          ic=ic+1;
          if is_even(ic)
           icv=icv+1;
           cceu(icv)=cced(ncr,ncc);
          end
         end
        end
        cce=cceu;
       end

       fiR=find(abs(cce)<=Rad_max & abs(cce)~=0);
       cces=cce;
       cce=cce(fiR);

%       figure, plot(cce,'r.')
%       grid, axis equal
%       fi=linspace(0,2*pi,200);
%       o=ones(size(fi))*Rad_max;
%       o1=ones(size(fi))*2*dcex;
%       hold on
%       plot(o.*exp(j*fi),'w')
%       plot(o1.*exp(j*fi),'y')
%       pausak

      case {'pch'}

       mx=ncex;
       my=ncey;
       ncex=(mx-1)*2+1;
       ncey=my;
       if is_even((mx-1)/2)==0
        Tyar=-1;
       else
        Tyar=1;
       end
       if Tyar==-1
        Tyar=2;
       end

       ncell=ncex*ncey;
       dcex=dxce/2;
       dcey=dyce*sqrt(3)/2;

       Rad_max=dcex*(mx-1);

       Rvetx=repmat(Rx,1,ncell);
       Rvety=repmat(Ry,1,ncell);
       Delta=repmat(Del,1,ncell);
       sh_type=repmat(sh_ty,1,ncell);

       fcex=([1:ncex]-ncex/2-.5)*dcex;
       fcey=([1:ncey]-ncey/2-.5)*dcey;
       Fcx=ones(size(fcey'))*fcex;
       Fcy=fcey'*ones(size(fcex));
       cced=Fcx+j*Fcy;

       if Tyar~=0
        ic=Tyar;
        icv=0;
        cceudu=ones(size(cced))*NaN;
        for ncc=1:ncex
         for ncr=1:ncey
          ic=ic+1;
          if is_even(ic)
           icv=icv+1;
           cceudu(ncr,ncc)=cced(ncr,ncc);
          end
         end
        end

        icdu=0;
        for ncr=1:ncey
         for ncc=1:ncex
          CEN=cceudu(ncr,ncc);
           if isnan(CEN)~=1
            icdu=icdu+1;
            if icdu~=1
             icv=icv+1;
             cceu(icv)=cceudu(ncr,ncc);
            end
           end
           if icdu==3
            icdu=0;
           end
         end
        end
        cce=cceu;

       end



       fiR=find(abs(cce)<Rad_max & abs(cce)~=0);
       cces=cce;
       cce=cce(fiR);

%       figure, plot(cce,'r.')
%       grid, axis equal
%       fi=linspace(0,2*pi,200);
%       o=ones(size(fi))*Rad_max;
%       o1=ones(size(fi))*2*dcex;
%       hold on
%       plot(o.*exp(j*fi),'w')
%       plot(o1.*exp(j*fi),'y')
%       pausak

      end

      end


      Pshd{1}=sh_type;
      Pshd{2}=Rvetx;
      Pshd{3}=Rvety;
      Pshd{4}=Delta;
      Pshd{5}=cce;
      Psh{icoso}=Pshd;

     end  %for kisor


    elseif sha==6
%       radii.array{pcounces,9}=lablar;
%       radii.array{pcounces,1}=labl;
%       radii.array{pcounces,2}=ori;
%       radii.array{pcounces,3}=shif;
%       radii.array{pcounces,4}=circle;
%       radii.array{pcounces,5}=D;
%       radii.array{pcounces,6}=d;
%       radii.array{pcounces,7}=Ry;
    isc=0;
    for isco=isor'
     isc=isc+1;
     PV{isc}.D=ary{isco}{5};
     PV{isc}.d=ary{isco}{6};
     PV{isc}.Ry=ary{isco}{7};
     PV{isc}.Rx=ary{isco}{8};
     PV{isc}.shape=ary{isco}{9};
     PV{isc}.orien=ary{isco}{2};
     PV{isc}.shif=ary{isco}{3};
     PV{isc}.Radd=ary{isco}{4};
     if length(ary{isco})>=10
      PV{isc}.grap=ary{isco}{10};
     else
      PV{isc}.grap=0;
     end
    end
     shape='sha_grat';
%     shape='sha_gjm';
%     disp(' mat_mix'), shape
%     keyboard
    end
   end

%disp(' mat_mix ver')
%     keyboard

    igintsa=igint;

%    if iLP==1
%     Kmat_LP0
%    else
%     Kmat_ve0
%    end

   shasav=sha;
   DOE=0;
   if sha==8
     sha=1
     DOE=1;
     fis=fis(1);
     adis=kcav0*radii.array{ish(fis),1}{10};
     Wdis=ones(size(adis));
     Wdis(2:2:end)=-1;
%     figure, plot(adis,Wdis,'wo')
     lxi=length(adis);
     aytot(:,8)=adis'/kcav0;
     ifii=-(1:lxi);
     ' solo un doe ! '
%' stop qui ', keyboard

   end

   if inume==0
%' stop qui kmat ', keyboard
    if iLP==1
     Kmat_LP
    else
    tic
     Kmat_ve
    toc
%    keyboard
    end
   else
%' stop qui kmat2 ', keyboard
    if iLP==1
     Kmat_LP2
    else
    tic
     Kmat_ve2
    toc
%    keyboard
    end
   end

%' stop qui dopo ', keyboard
   if shasav==8
     sha=8
     DOE=0;
%' stop qui dopo ', keyboard
   end


    if (sha==6 | sha==5 | sha==3) & ifp~=-4 & ipolar~=0
     map(reshape(Kosp(:,:,1),si2)),
     map(reshape(Kosm(:,:,1),si2)),
     DK=Kosm-Kosp;
     map(reshape(DK(:,:,1),si2)),
     drawnow
    end
%    keyboard
    igint=igintsa;

    if ireturn==0
     sK=length(size(Kosp));

     MKosp{shad}=Kosp;
     MKosm{shad}=Kosm;
%     disp('Mkos'), keyboard
     if iztm==1
      MKoszp{shad}=Koszp;
      MKoszm{shad}=Koszm;
     end
    else
      grat_app
    end  %ireturn


%  end

 end  % shai==sha
 disp(' loop sha: '), sha
% clear KAp KAm KAzp KAzm
 clear Kosp Kosm Koszp Koszm
end

%' MKosp ', keyboard
% rid_bas1
% rid_bas3
% rid_bas2




if length(iredmat)==0
 iredmat=0;
end



%fir=find(shavet==6);
%fir=find(shavet==6 & dv~=0);
fir=find(shavet==6 | abs(shavet)==8);
if ikr==0
 fir=[];
end
%' fir=find(shavet==6); ', keyboard


if nucomp==0
 fir=[];
end


%if inume==1
%
% if ifr==1
%  if length(fir)==0
%   rid_uu
%  else
%   rid_reta
%  end
% end
%
%else
%
% if ifr==1
%  if length(fir)==0
%   rid_uu
%  end
% end
%
%end


 if ifr==1
  if length(fir)==0
   rid_uu
  end
 end
 if ifp==-10
  ' rid_uu ', keyboard
 end

fish=find(shavet~=0 & dvc~=0);
[shai,ishd]=sort(shavet(fish));
fisha=find(diff([0; shai])~=0);
shailoop=shai(fisha)';
shailoop=abs(shailoop);

for sha=shailoop
%  if sha==7
%   sha=1;
%  end
  Kplot1=MKosp{sha};
  Kplot2=MKosm{sha};
  if inume==0 & iredmat==1
   Kplot1=Kplot1(Pusa,Pusa,:);
   Kplot2=Kplot2(Pusa,Pusa,:);
%   MKosp{sha}=Kplot1;
%   MKosm{sha}=Kplot2;
%   if iztm==1
%    Kploz1=MKoszp{sha};
%    Kploz2=MKoszm{sha};
%    Kploz1=Kploz1(Pusa,Pusa,:);
%    Kploz2=Kploz2(Pusa,Pusa,:);
%    MKoszp{sha}=Kploz1;
%    MKoszm{sha}=Kploz2;
%    clear Kploz1 Kploz2
%   end
  end
%  ' prima di anmat ', keyboard
   if ianmat==1
    an_mat
%  ' dopo anmat ', keyboard
   else
      iaccp =[];
      iaccm =[];
      PUridm=[];
      PUridp=[];
   end
%   ' sha ', pausak
      Iaccp{sha}=iaccp;
      Iaccm{sha}=iaccm;
      PUriM{sha}=PUridm;
      PUriP{sha}=PUridp;
end


%' dopo rid ext nuovo', keyboard


% calcolo i coeff. di accop. termici
KApsa=KAp;
KAmsa=KAm;
%'% calcolo i coeff. di accop. termici', keyboard, keyboard
%
    if (igau==4 & length(yiT)>1)
     disp(' calcolo i coeff. di accop. termici ')
%     keyboard

     sha=1;
     adis=xiT*kcav;
     lxi=length(adis);
     igainshape=1;

     if iLP==1
      Kmat_LP
     else
      Kmat_ve
     end

%' dopo Kmat ', keyboard
% lo salto perche` e` gia` fatto dentro a *_ci.m
ilosal=1;
if ilosal==0

     s=size(yiT);
     if iLP==1

      for iz=1:s(2)
       Ktoiie=0;
       indis=0;
       for ndis=1:lxi
        indis=indis+1;
        Kdum=reshape(Kti(ndis,:,:),si);
        Ktoiie=Ktoiie+Kdum*yiT(indis,iz);
       end  %ndis
       if isi==0
        KAp=diag(ZEv.*KKv)*(Ktoiie);
       else
        diae=diag(sqrt(ZEv.*KKv));
        KAp=diae*(Ktoiie)*diae;
       end
       KAp=KAp(Pusa,Pusa);
       KTem_p(:,:,iz)=KAp;
       KTem_m(:,:,iz)=KAp;

      end
%' dopo Kmat LP', keyboard

     else  %iLP


      for iz=1:s(2)
       Ktoiie=0;
       Ktoije=0;
       Ktoiim=0;
       Ktoijm=0;
       Ktoiiz=0;
       indis=0;
       for ndis=1:lxi
        indis=indis+1;
        Kdum=reshape(Ktoiiei(ndis,:,:),si);
        Ktoiie=Ktoiie+Kdum*yiT(indis,iz);
        Kdum=reshape(Ktoijei(ndis,:,:),si);
        Ktoije=Ktoije+Kdum*yiT(indis,iz);
        Kdum=reshape(Ktoiimi(ndis,:,:),si);
        Ktoiim=Ktoiim+Kdum*yiT(indis,iz);
        Kdum=reshape(Ktoijmi(ndis,:,:),si);
        Ktoijm=Ktoijm+Kdum*yiT(indis,iz);
        if iztm==1
         Kdum=reshape(Ktoiizi(ndis,:,:),si);
         Ktoiiz=Ktoiiz+Kdum*yiT(indis,iz);
        end
       end  %ndis
       if isi==0
        KEEp=diag(ZEv.*KKv)*(Ktoiie);
        KEMp=diag(ZEv.*KKv)*(Ktoije);
        KMEp=diag(ZMv.*KKv)*(Ktoijm);
        KMMp=diag(ZMv.*KKv)*(Ktoiim);
       else
        diae=diag(sqrt(ZEv.*KKv));
        diam=diag(sqrt(ZMv.*KKv));
        diae(1,1)=diae(1,1)*fatde1;
        diam(1,1)=diam(1,1)*fatde1;
        KEEp=diae*(Ktoiie)*diae;
        KEMp=diae*(Ktoije)*diam;
        KMEp=diam*(Ktoijm)*diae;
        KMMp=diam*(Ktoiim)*diam;
       end
       KAp=[KEEp KEMp; KMEp KMMp];
       KAm=[KEEp KEMp; KMEp KMMp];
       if iztm==1
        if isi==0
         KMMpz=diag(0.5*ZMv.*KKv.*KKv./bev)*(Ktoiiz)*diag(KKv./bev);
        else
         diam=sqrt(ZMv.*KKv);
         diam(1)=diam(1)*fatde1;
         KMMpz=diag(0.5*diam.*KKv./bev)*(Ktoiiz)*diag(KKv./bev.*diam);
        end
        KAzp=[Kzer Kzer; Kzer KMMpz];
       else
        KAzp=0;
       end
       KAp=KAp(Pusa,Pusa);
       KAm=KAm(Pusa,Pusa);
       KAzp=KAzp(Pusa,Pusa);
       KTem_p(:,:,iz)=KAp;
       KTem_m(:,:,iz)=KAm;
       if iztm==1
        KTem_z(:,:,iz)=KAzp;
       end
      end
%' dopo Kmat Vett', keyboard

     end
else
 si2=[length(Pusa) length(Pusa)];
 if iLP==1
  KTem_p=KTemp(Pusa,Pusa,:);
  KTem_m=KTemp(Pusa,Pusa,:);
 else
  KTem_p=KTemp(Pusa,Pusa,:);
  KTem_m=KTemp(Pusa,Pusa,:);
  KTem_z=KTempz(Pusa,Pusa,:);
 end

end  %ilosal
%' ilosal ' , keyboard
    else

       KTem_p=0;
       KTem_m=0;
       KTem_z=0;

    end

KAp=KApsa;
KAm=KAmsa;
%'quia', keyboard

% calcolo i coeff. di accop. della parte attiva
%

    if ( igau==0 | (igau==4 & length(yiN)<=1) )

     if ~exist('KAp')
      disp(' errore 1 in mat_mix: manca KA area costante! ')
      pausak
     end

    elseif (abs(igau)==4 & length(yiN)>1)

disp(' calcolo i coeff. di accop. della parte attiva ')
%keyboard
%     sha=1;
%     adis=xiN*kcav0;
%     lxi=length(adis);
%     igainshape=1;

     fiat=find(iauto(:,1)==2);
     sha=shavet(fiat,1);
     shape=sha_fun(sha);
     rassi=radii.b(fiat,1)/radii.a(fiat,1);
     adis=xiN*kcav;
     bdis=xiN*kcav*rassi;
%     aloc=adis;
%     bloc=bdis;
%     adis=xiN*kcav0*rassi;
%     bdis=xiN*kcav0;
     pdi=ones(size(xiN))*radii.c(fiat,1);
     lxi=length(adis);
     igainshape=1;
% 'Kmat attivo', keyboard

     if iLP==1
%      'Kmat_LP', keyboard

      Kmat_LP
     else
%      'Kmat_ve', keyboard
      Kmat_ve
     end


disp(' dopo calcolo i coeff. di accop. della parte attiva ')
%keyboard
    if sha==1  %altri casi gia fatti in Kmat
     if iLP==1

       Ktoiie=0;
       indis=0;
       for ndis=1:lxi
        indis=indis+1;
        Kdum=reshape(Kti(ndis,:,:),si);
        Ktoiie=Ktoiie+Kdum*yiN(indis);
       end  %ndis


       if ipolar~=0
        if isi==0
         KAp=diag(ZEv.*KKv)*(Ktoiie);
        else
         diae=diag(sqrt(ZEv.*KKv));
         KAp=diae*(Ktoiie)*diae;
        end
       else
        KAp=Ktoiie;
       end

       if ianti_gui==1 & exist('yiN_ag')
        Kag=0;
        indis=0;
        for ndis=1:lxi
         indis=indis+1;
         Kdum=reshape(Kti(ndis,:,:),si);
         Kag=Kag+Kdum*yiN_ag(indis);
        end  %ndis
        Ktoiiea=Kag;
        if ipolar~=0
         if isi==0
          KA_ag=diag(ZEv.*KKv)*(Ktoiiea);
         else
          diae=diag(sqrt(ZEv.*KKv));
          KA_ag=diae*(Ktoiiea)*diae;
         end
        else
         KA_ag=Ktoiie;
        end
        KA_ag=KA_ag(Pusas,Pusas);
       end

       KAp=KAp(Pusas,Pusas);
       KAm=KAp;
%       'qui  LP', keyboard

     else  %iLP


       Ktoiie=0;
       Ktoije=0;
       Ktoiim=0;
       Ktoijm=0;
       Ktoiiz=0;
       indis=0;
       for ndis=1:lxi
        indis=indis+1;
        Kdum=reshape(Ktoiiei(ndis,:,:),si);
        Ktoiie=Ktoiie+Kdum*yiN(indis);
        Kdum=reshape(Ktoijei(ndis,:,:),si);
        Ktoije=Ktoije+Kdum*yiN(indis);
        Kdum=reshape(Ktoiimi(ndis,:,:),si);
        Ktoiim=Ktoiim+Kdum*yiN(indis);
        Kdum=reshape(Ktoijmi(ndis,:,:),si);
        Ktoijm=Ktoijm+Kdum*yiN(indis);
        if iztm==1
         Kdum=reshape(Ktoiizi(ndis,:,:),si);
         Ktoiiz=Ktoiiz+Kdum*yiN(indis);
        end
       end  %ndis

       if isi==0
        KEEp=diag(ZEv.*KKv)*(Ktoiie);
        KEMp=diag(ZEv.*KKv)*(Ktoije);
        KMEp=diag(ZMv.*KKv)*(Ktoijm);
        KMMp=diag(ZMv.*KKv)*(Ktoiim);
       else
        diae=diag(sqrt(ZEv.*KKv));
        diam=diag(sqrt(ZMv.*KKv));
        diae(1,1)=diae(1,1)*fatde1;
        diam(1,1)=diam(1,1)*fatde1;
        KEEp=diae*(Ktoiie)*diae;
        KEMp=diae*(Ktoije)*diam;
        KMEp=diam*(Ktoijm)*diae;
        KMMp=diam*(Ktoiim)*diam;
       end
       KAp=[KEEp KEMp; KMEp KMMp];
       KAm=[KEEp KEMp; KMEp KMMp];
       KAp=KAp(Pusas,Pusas);
       KAm=KAm(Pusas,Pusas);

       if iztm==1
        if isi==0
         KMMpz=diag(0.5*ZMv.*KKv.*KKv./bev)*(Ktoiiz)*diag(KKv./bev);
        else
         diam=sqrt(ZMv.*KKv);
         diam(1)=diam(1)*fatde1;
         KMMpz=diag(0.5*diam.*KKv./bev)*(Ktoiiz)*diag(KKv./bev.*diam);
        end
        KAzp=[Kzer Kzer; Kzer KMMpz];
        KAzp=KAzp(Pusas,Pusas);
        KAzm=KAzp;
       else
        KAzp=0;
        KAzm=0;
       end

%       'qui  VE prima', keyboard

       if ianti_gui==1 & exist('yiN_ag')

         Ktoiie=0;
         Ktoije=0;
         Ktoiim=0;
         Ktoijm=0;
         Ktoiiz=0;
         indis=0;
         for ndis=1:lxi
          indis=indis+1;
          Kdum=reshape(Ktoiiei(ndis,:,:),si);
          Ktoiie=Ktoiie+Kdum*yiN_ag(indis);
          Kdum=reshape(Ktoijei(ndis,:,:),si);
          Ktoije=Ktoije+Kdum*yiN_ag(indis);
          Kdum=reshape(Ktoiimi(ndis,:,:),si);
          Ktoiim=Ktoiim+Kdum*yiN_ag(indis);
          Kdum=reshape(Ktoijmi(ndis,:,:),si);
          Ktoijm=Ktoijm+Kdum*yiN_ag(indis);
          if iztm==1
           Kdum=reshape(Ktoiizi(ndis,:,:),si);
           Ktoiiz=Ktoiiz+Kdum*yiN_ag(indis);
          end
         end  %ndis

         if isi==0
          KEEp=diag(ZEv.*KKv)*(Ktoiie);
          KEMp=diag(ZEv.*KKv)*(Ktoije);
          KMEp=diag(ZMv.*KKv)*(Ktoijm);
          KMMp=diag(ZMv.*KKv)*(Ktoiim);
         else
          diae=diag(sqrt(ZEv.*KKv));
          diam=diag(sqrt(ZMv.*KKv));
          diae(1,1)=diae(1,1)*fatde1;
          diam(1,1)=diam(1,1)*fatde1;
          KEEp=diae*(Ktoiie)*diae;
          KEMp=diae*(Ktoije)*diam;
          KMEp=diam*(Ktoijm)*diae;
          KMMp=diam*(Ktoiim)*diam;
         end
         KAp_ag=[KEEp KEMp; KMEp KMMp];
         KAm_ag=[KEEp KEMp; KMEp KMMp];
         KAp_ag=KAp_ag(Pusas,Pusas);
         KAm_ag=KAm_ag(Pusas,Pusas);

         if iztm==1
          if isi==0
           KMMpz=diag(0.5*ZMv.*KKv.*KKv./bev)*(Ktoiiz)*diag(KKv./bev);
          else
           diam=sqrt(ZMv.*KKv);
           diam(1)=diam(1)*fatde1;
           KMMpz=diag(0.5*diam.*KKv./bev)*(Ktoiiz)*diag(KKv./bev.*diam);
          end
          KAzp_ag=[Kzer Kzer; Kzer KMMpz];
         else
          KAzp_ag=0;
         end
         KAzp_ag=KAzp_ag(Pusas,Pusas);
         KAzm_ag=KAzp_ag;


%        Kag=0;
%        indis=0;
%        for ndis=1:lxi
%         indis=indis+1;
%         Kdum=reshape(Kti(ndis,:,:),si);
%         Kag=Kag+Kdum*yiN_ag(indis);
%        end  %ndis
%        Ktoiiea=Kag;
%        if ipolar~=0
%         if isi==0
%          KA_ag=diag(ZEv.*KKv)*(Ktoiiea);
%         else
%          diae=diag(sqrt(ZEv.*KKv));
%          KA_ag=diae*(Ktoiiea)*diae;
%         end
%        else
%         KA_ag=Ktoiie;
%        end

       end

     end
    else  %sha
%    ' verifica prima', keyboard

        if iLP==0
         KAp_ag=KAp_ag(Pusas,Pusas);
         KAm_ag=KAm_ag(Pusas,Pusas);
         KAp=KAp(Pusas,Pusas);
         KAm=KAm(Pusas,Pusas);
         if iztm==1
          KAzp_ag=KAzp_ag(Pusas,Pusas);
          KAzm_ag=KAzm_ag(Pusas,Pusas);
          KAzp=KAzp(Pusas,Pusas);
          KAzm=KAzm(Pusas,Pusas);
         end
        else
         KAp_ag=KAp_ag(Pusas,Pusas);
         KAm_ag=KAm_ag(Pusas,Pusas);
         KAp=KAp(Pusas,Pusas);
         KAm=KAm(Pusas,Pusas);
        end

    end  %sha
%    ' verifica ', keyboard
%      if ifp==-10
%       'qui fine K attivo', keyboard
%      end

    else
     disp(' errore 2 in mat_mix: manca KA! ')
     keyboard
    end

     if ifp~=-4 & ifr==1
      map(KAp)
      title(' KAp ')
      pausak
     end

% calcolo i coeff. di accop. dell'anisotropia e Idelta
%
%sha
%keyboard
 Kmat_any


%Kd=MKosp{6};
%' iredmat ', keyboard

if iredmat==1
% Idelta=Idelta(Pusas,Pusas);
% if iztm==1
%  Ideltaz=Ideltaz(Pusas,Pusas);
% end
%' sav', keyboard
 KAp=KAp(Pusasav,Pusasav);
 KAm=KAm(Pusasav,Pusasav);
 if iztm==1
  KAzp=KAzp(Pusasav,Pusasav);
  KAzm=KAzm(Pusasav,Pusasav);
 end
 sM=size(MKosp);
 MKospd=MKosp;
 MKosmd=MKosm;
 if iztm==1
  MKoszpd=MKoszp;
  MKoszmd=MKoszm;
 end
 clear MKosp MKosm   MKoszp MKoszm
%' cosa succede?? ', keyboard
 for ki=1:sM(2)
  sMi=size(MKospd{ki});
  if sum(sMi)>0
   if length(sMi)==3
    kjM=sMi(3);
   else
    kjM=1;
   end
%   MKosp{ki}=[];
%   MKosm{ki}=[];
%   MKoszp{ki}=[];
%   MKoszm{ki}=[];
   for kj=1:kjM
    MKosp{ki}(Pus,Pus,kj)=MKospd{ki}(Pusasav,Pusasav,kj);
    MKosm{ki}(Pus,Pus,kj)=MKosmd{ki}(Pusasav,Pusasav,kj);

%    if ianmat==0 & ki==1
%   ' ianmat', keyboard
%     an_mat
%      Iaccp{1}=iaccp;
%      Iaccm{1}=iaccm;
%      PUriM{1}=PUridm;
%      PUriP{1}=PUridp;
%    end
     if iztm==1
      MKoszp{ki}(Pus,Pus,kj)=MKoszpd{ki}(Pusasav,Pusasav,kj);
      MKoszm{ki}(Pus,Pus,kj)=MKoszmd{ki}(Pusasav,Pusasav,kj);
     end
   end
  end
 end
end


Kosm=MKosm;
Kosp=MKosp;
if iztm==1
 Koszm=MKoszm;
 Koszp=MKoszp;
end


if ifp~=-4
 disp(' fine mat_mix '),
end
%keyboard
%keyboard

if ifp>-4
 keyboard
end

clear KA KAmdu KAmsub KApdu KApsub KAz KAz1 KAzmdu KAzmsub KAzmz
clear KAzpdu KAzpsub KBz1 KEEp KEMp KMEp KMMp KMMpz Ktoa Ktoa0
clear KTemp KTempz Kadia Kany Kany1 Kany2 Kanyd Kanyd1
clear Kanyd2 Kd1e Kd1m Kde Kdia Kdm Kdz Kiie Kiim Kiiz Kije Kijm Kos  Kosmsub
clear Kospsub Kosz  Koszmsub  Koszpsub Kpm Kpp Ksm Ksp
clear Ktoa0v Ktoa1 Ktoad Ktoiie Ktoiiei Ktoiim Ktoiimi Ktoiiz Ktoiizi Ktoije
clear Ktoijmi Kzer Ktoijei Ktoijm Kset Vset
clear B A
clear Kplot1 Kplot2 Kp1 Kpp1 Kp2 Kpp2
clear KMMmz KApsa KAmsa

clear Fi1 Fi2 Fi3 Fi4 Fi5 Fi6 Ide
clear MKosp MKosm MKoszp MKoszm


[Tfas,Ffas]=eltime(Tfas,Ffas,-10);



if ifp==-10
 disp(' fine mat_mix dopo clear '),
% keyboard
end

% disp(' fine mat_mix dopo clear '),
% keyboard

%%% inizializzazione variabili
pack

P=ones(2*lKA,2*lKA)*(1+i);

if length(inuo_bas)>0
%' qui vedo ', keyboard
 global  Oo Odu
%' qui vedo '
 Odu=P;
end


Oo=ones(size(IdeOon))*(1+i);
Tc=ones(size(IdeOon));
T=Tc;
Tcm=Tc;
Tdu=Tc;
Tdum=Tc;
Tmirro=Tc;
Tco=Tc;
%Tstof=Tc;
Mn=Tc;
Mi=Tc;
Moi=Tc;
Mon=Tc;
KOz=ones(lKA,lKA)*(1+i);
KOt=KOz;
Kr=KOz;
Krz=KOz;
Kost=KOz;
Kostz=KOz;
Kr=KOz;
Krz=KOz;
Ksu=KOz;
pI=KOz;
pIz=KOz;

' fine  set variabili ',
mem_used
%keyboard
if exist('si2sav')
  si2=si2sav;
end


%' dopo mat', keyboard
%
%if length(inuo_bas)>0
%
% inuo_bas
%'nuo_bas', keyboard
%
% nuo_bas
%end

%global xprod yprod
%xprod=Tc;
%yprod=Tc;
