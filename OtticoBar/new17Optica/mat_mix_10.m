%emme='emme_ult';
%load matmix
%save matmix
'mat_mix 10',  keyboard
 KANzp=0;
 KANzm=0;
icomp_anymat=0;
Madd=0;
shaC=0;
segem=-1;   %mette i segni come Love. =1 usa i segni TM cambiati, come fatto dall'inizio
%segem=1;   %mette i segni come Love. =1 usa i segni TM cambiati, come fatto dall'inizio

Madd0=0;
iclo=0;
igainshapeT=0;
anyret=zeros(size(dv));
pros=anyret;
nv=nv_sa;
for kseg=1:length(nv)
duse=find(real(nv(kseg,:)<0));
 if length(duse)>0
  pros(kseg,1)=-duse;
  if duse>1 & duse<length(find(nv(kseg,:)~=0))
   pros(kseg,2)=-duse-1;
  end
 else
  pros(kseg,1)=1;
 end
end 
fiseg=find(pros<0);
%shavet(fiseg)=-shavet(fiseg);
%'trovo raggi anisotropi', keyboard
if length(fiseg)>0             %specchio curvo con reticolo
 rad_du=rad(fiseg);
 raso=sort(rad_du);
 fiz=find(diff([raso; 0])~=0);
 rad_any=raso(fiz);
 if rad_any(1)>0 & igr_app>0
  rad_any=[0; rad_any];
 end
 
 PAlei=PAle{1};
 P=PAlei;
 Psalva=PAlei;
 PAnri=PAnr{1};
 nlen=PAnri;
 pos=PAlei.pos;
 UD=PAlei.UD;
 

 ngra=ngvero;

 t1=P.d;
 t2=P.LA-P.d;
 tt=P.LA;
 rad=radii.a;
 
 if igr_app>0
  er1=ngra(1)^2;
  er2=ngra(2)^2;
  n_ve=tt*er1*er2/(t2*er1+t1*er2);
  n_pa=(t1*er1+t2*er2)/tt;
  n_me=((n_ve+n_pa)/2);
  n_di=((n_pa-n_ve)/2);
  if P.orien==0
   ey=n_ve;
   ex=n_pa;
  else
   ex=n_ve;
   ey=n_pa;
  end
  Delta_eps=(ex-ey)/2; 
  
  anyret(fiseg)=Delta_eps;
  
%  'anyret', keyboard
  fimin=find(real(nv)<0);
  nv(fimin)=sqrt(n_me);
  icomp_anymat=1;
 else  %igr_app
  fimin=find(real(nv)<0);
  nv(fimin)=ngra(1);
  icomp_anymat=2; 
  fi0=find(rad_any==0);
  RAmax=PAlei.D;
  rad_any(fi0)=RAmax;
  ragra=rad(fiseg);
  fif0=find(ragra==0);
   if length(fif0)>0
    fifset=fiseg(fif0);
    rad(fifset)=RAmax;
   end 
   aitn=rad(put,:).';
   aib=rad(pub,:).';
   for kse=fiseg'
    rai=rad(kse,:);
    if rai==RAmax
      nvii=nv(kse,:);
      fizero=find(real(nvii)~=0);
      fisos=fizero+1;
      nv(kse,fisos)=1;
    end
   end
%  rad_any=sort(rad_any);
 end
 if ifp==-10
  'qui aniret in matmix', keyboard
  'qui aniret in matmix', keyboard
 end
else
 KANp=0;
 KANm=0;
end  %fiseg


%if exist('anyr')==0
 anyr.t=anyret(put,:);
 anyr.b=anyret(pub,:);
 nitn=nv(put,:).';
 nib =nv(pub,:).'; 
 prosn=pros(put,:).';
 prosb =pros(pub,:).';  
%end  

if ifp==-10
%'qui mat_mix ingresso', keyboard
end

%global inuo_bas iredmat nucomp
global inuo_bas nucomp
clear Kosp Kosm Koszp Koszm
if length(pMu0u)>0
%clear pMu0u
%clear global pMu0u
end
%' entro mix', keyboard


%if length(inuo_bas)==0
% inuo_bas=0;
%end
inuo_bas=1;
%inuo_bas=0;
global irid_bas
global igr_app
%if igr_app==1
if igr_app>=1
 inuo_bas=1;
 irid_bas=1;
% irid_bas=0;
 lKAn=length(KK)*2;
 IdeOon=diag(ones(2*lKAn,1));
 Pus=1:lKAn;
 Pust=1:lKAn*2;
end

ianmat=1;
ianmat=0;
if pasnu==1
 ianmat=0;
end
global igr_app
if igr_app>=1
%if igr_app==1
 ianmat=0;
end
%'ian', keyboard
global iredmat
if inuo_bas==1
  iredmat=1;
else
 iredmat=0;
 nucomp=0;
end
%iredmat=0;
if ifr==1
if iLP==1
 Pusas=1:(numodiacc+1)*length(KK);
else
 Pusas=1:(numodiacc+1)*length(KK)*2;
end
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
fis0=find(radii.a==0 & shavet~=6);
sizsh=size(shavet);
dvc=repmat(dv,1,sizsh(2));
shavet(fis0)=0;

shavdu=reshape(shavet,prod(size(shavet)),1);
dvcv=reshape(dvc,prod(size(shavet)),1);
[shad,ipush]=sort(shavdu);
fino0=find(shad~=0 & dvcv(ipush)~=0);
ish=ipush(fino0);
shai=shavet(ish);
%' primo shai ', keyboard



ax=radii.b(ish);
ay=radii.a(ish);

pad=radii.c(ish);
sha_ish=shavet(ish);
aryd=radii.array;
ary=repmat(aryd,1,size(shavet,2));
ary_s=ary(ish);
ary=ary_s;


%firet=find(shavet(ish)==6 & dv();
firet=find(shavet(ish)==6);
%'firet ', keyboard

if length(firet)>=1
 for kfir=firet'
  icel=ish(kfir);
  DU=radii.array{icel};
  pad(kfir)=DU{1};

 end

%end  % prima era qui


 ico=0;

 for duish=1:length(aryd)
  ico=ico+1;
  du=aryd{ico}{1};
  if length(du)==0
   du=0;
  end
  ty(ico,:)=du*ones(1,size(shavet,2));
 end

end
%end  % prima era qui


 ity=0;
if length(find(shavet>4 & shavet~=7 & shavet~=6))>0
%if length(find(abs(shavet)>4 & shavet~=7 & shavet~=6))>0
 ishve=find(shavet>4 & shavet~=7);
 ishve=find(shavet(:,1)>0 & shavet(:,1)~=7);
% ishve=find(abs(shavet)>0 & shavet~=7);
% ary=radii.array(ishdu);
% aryd=radii.array(ishve);
% 'passso invece qui', keyboard
% aryd=radii.array(ishdu_a);
 ico=0;
% for duish=ishdu'
 for duish=ishve'
  ico=ico+1;
  ty(ico,1)=ary{ico}{1}*ones(1,size(shavet,2));
 end
 ity=1;
end

% aryd=radii.array(ishdu_a);
% ary=repmat(aryd,1,size(shavet,2));

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
if length(fir)==0
 if igr_app>0
  fir=0;
 end
end
if icomp_anymat==2
 fir=0;
end
%fir
%keyboard
%' fir=find(shavet==6); ', keyboard
%' mat_mix ', pausak


%' prima riretu', keyboard
if length(nucomp)==0
 nucomp=0;
end
if nucomp==0
 inume=0;
else
 inume=1;
end
if fir~=0
 inume=1;
end

if inume==1
 if ifr==1 & iredmat==1 
% if length(fir)>0
  ri_retu
  si2sav=si2;
  if ifp~=-4
% figure, plot(KKv(Pusas),'.')
% Pus0=Pusas;
    figure, plot(KKv(Pus0),'.')
    ' qui iNM', keyboard
    if iclo==1
     close
    end
  end
% end
 end
%kmat_any
%' kmat_any', keyboard

end
% if ifp~=-4 & exist('Pus0')
 if ifp==-40 & exist('Pus0')
   Pud=Pus0;
  figure,
  subplot(211)
  plot(KKv(Pud),real(Gas(Pud)),'.',KKv(Pud),real(Gad(Pud)),'.'),
  subplot(212)
  plot(KKv(Pud),imag(Gas(Pud)),'.',KKv(Pud),imag(Gad(Pud)),'.'),
  pausak

   if iclo==1
    close
   end

 end

ifipu=ifil(ish);
ifipud=ifipu;


fisha=find(diff([0; shai])~=0);
shailoop=shai(fisha)';

shailoop=abs(shailoop);
shailoopu=abs(shailoop);

shai=abs(shai);

%shailoop=4
%' sahiloop', keyboard
for shad=shailoop
ireturn=0;
 if shad<0
  shad=-shad;
 end
 sha=abs(shad);
 if ifp~=-4
 disp([' mat_mix: sha = ',num2str(sha)])
 end
% ' saha ', keyboard
 fis=find(shai==sha);
   if sha==7 | sha==4
    if sha==7
     sha=1;
    end
%    ' sha = 7 in mat_mix', keyboard
%    if length(find(Raxv~=1))>0
%     shape='sha_rhom';
%     sha=4;
%    end
    if length(find(shavet==7))>0
     shaC=7;
    else 
     shaC=0;
    end
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

  axi=ax(isor);
  ayi=ay(isor);
  pdi=pad(isor);
%disp(' mat_mix interno'), keyboard
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

'cce', keyboard
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
     if length(ary{isco})>=11
      PV{isc}.Rext=ary{isco}{11};
     else
      PV{isc}.Rext=0;
     end
    end
     shape='sha_grat';
%     shape='sha_gjm';
%     disp(' mat_mix'), shape
%     keyboard
    end
   end

%disp(' mat_mix ver'),     keyboard

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
%     axtot(:,8)=0;
     ifii=-(1:lxi);
     ' solo un doe ! '
   end
%    ' stop qui kmat_ve ', keyboard
   if inume==0
  if ifp==-10
%    ' stop qui kmat_ve ', keyboard
   end
    if iLP==1
     Kmat_LP
    else
    tic
     Kmat_ve
    toc
%    keyboard
    end
   else
   if ifp==-10
%    ' stop qui kmat_ve2 ', keyboard
   end
    if iLP==1
     Kmat_LP2
    else
    tic
     Kmat_ve2
    toc
%    keyboard
    end
   end

   if shasav==8
     sha=8
     DOE=0;
   end


    if (sha==6 | sha==5 | sha==3) & ifp~=-4 & ipolar~=0
      ' >>>>>>>>reshape'
      keyboard
      keyboard

     map(Kosp(:,:,1)),
     map(Kosm(:,:,1)),
     DK=Kosm-Kosp;
     map(DK(:,:,1)),
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
% disp(' loop sha: '), sha
% keyboard
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
if icomp_anymat==2
 fir=0;
end
%if length(fir)==0
% if igr_app>0
%  fir=0;
% end
%end

%fir
%keyboard
%' fir=find(shavet==6); ', keyboard


if nucomp==0 & igr_app==0
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
% end7
%
%end


%' prima rid_uu ', keyboard
 if ifr==1
  if length(fir)==0 
   rid_uu
  end
 end
%' rid_uu ', keyboard

%fish=find(shavet~=0 & dvc~=0);
%[shai,ishd]=sort(shavet(fish));
%fisha=find(diff([0; shai])~=0);
%shailoop=shai(fisha)';
%shailoop=abs(shailoop);
%'shailoopu', keyboard
for sha=shailoopu
%  if sha==7
%   sha=1;
%  end
  Kplot1=MKosp{sha};
  Kplot2=MKosm{sha};
  if ~exist('Pusa')
   Pusa=Pusas;
  end

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
  ' prima di anmat ', keyboard
   if ianmat==1
    an_mat
%    an_mat_old
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
if iplan==0
KApsa=KAp;
KAmsa=KAm;
end
%'% calcolo i coeff. di accop. termici', keyboard, keyboard

if ifp==-10
 'fine Mat_mix'
 keyboard
end
%
    if (igau>=4 & length(yiT)>1)
     disp(' calcolo i coeff. di accop. termici ')
%     keyboard

     tic
     sha=1;
%     sT=ismat(T);

     if sT==2
      adis=xiT*kcav;
      if length(find(shavet==4))>1
       ips=find(shavet==4 & real(nv(:,2))<2 & abs(imag(nv(:,2)))<.01 ...
               & real(nv(:,2))>1);
       r_dr=radii.b(ips)/radii.a(ips); %attenzione: a b = y x

       icir=1;
       if icir==0
        sha=4;
        pdi=pdi(1)*ones(size(bdis));
       else
        r_dr=1;
        sha=1;
         pdi=1e-4*ones(size(bdis));
       end

%       ' qui per set parametri ', keyboard
        iresK=0;
%       bdis=(xiT+dr)*kcav;
        bdis=adis*r_dr;
       end
     elseif sT==3  %sT
        sha=4;
        iresK=1;
        global Pterm PcTerm
        zcut=PcTerm.z;
        zz=PcTerm.zv;
        roT=PcTerm.xv;
        Nz=length(Pterm);
        sP=size(Pterm{1});
        Nr=sP(2);
        xd=[]; yd=[]; pd=[]; vd=[]; Tdis=[];
        puP0=2:Nr;
        puP=1:Nr-1;
        for kz=1:Nz

         Pm=(Pterm{kz}(:,puP)+Pterm{kz}(:,puP0))/2;
%         xd=[xd; Pterm{kz}(1,puP)'];
%         yd=[yd; Pterm{kz}(2,puP)'];
%         pd=[pd; Pterm{kz}(3,puP)'];
%         vd=[vd Pterm{kz}(4,:)'];
         xd=[xd; Pm(1,:)'];
         yd=[yd; Pm(2,:)'];
         pd=[pd; Pm(3,:)'];
         vd=[vd Pterm{kz}(4,:)'];
         Tdis=[Tdis; fliplr(Pterm{kz}(4,:))];
        end
        Tdis=Tdis';
        xiT=xd+yd+pd;
        bdis=xiT*kcav;
        yiT=xd-yd+pd;
        adis=yiT*kcav;
        pdi=pd;

%       ' spnmo in ,at_mix', keyboard
%        ro_inf=(max(yiT)+3);
%        [du,fibo]=min(ro_inf-roT);

%         fiend=
        Tend=T(end,:,1);
        Tdu=spline(zz,Tend,zdis);
        Tdis(end,:)=Tdu;
        vd(1,:)=Tdu;
        dp=diff(vd);
        yiT=dp;
        if ifp~=-4
         figure, plot(cumsum(yiT',2),'.')
         ' somo in mat_mix', pausak
        end
     end


       if ifp~=-4 & is_even(mm+iLP)==0
        ' verifica discretizzazione '
        puzd=[1:length(zedis)];
        puz=[];
        for k=puzd
         [du,pp]=min(abs(zedis(k)-zeta));
         puz=[puz pp];
        end
%        zes=sum(yiT(:,puzd))'*ones(1,length(yiT));
        zes=sum(yiT(:,puzd))'*ones(1,size(yiT,1));
        if sT==2
         xiP=xiT;
         yiP=zes-cumsum(yiT(:,puzd),1)';
        else
         xiP=reshape(xiT',size(yiT));
%         xiP=xiP(:,puzd)';
%         yiP=cumsum(yiT(:,puzd),1)';
         xiP=xiP(:,puzd);
         yiP=cumsum(yiT(:,puzd),1);
        end
        sdi=size(Tdis);
        yad=(ones(sdi(1)-1,1)*Tdis(end,:));
        if sT==2
         yad=yad';
        end
%        figure,
%        plot(ro_inT,T(1:end-1,puz)-ones(length(T)-1,1)*T(end-1,puz))
%        hold on
%        plot(xiP,yiP,'.')
%        pausak
        figure,
        global Temper_sav
        if length(Temper_sav)>0
         plot(ro_inT,Temper_sav(1:length(ro_inT),puz))
        else
         plot(ro_inT,T(1:length(ro_inT),puz))
        end
        hold on
        plot(xiP,yiP+yad,'.')
        pausak
       end

     lxi=length(adis);
     igainshape=1;
     igainshapeT=1;
     ifp_salvo=ifp;
%     ifp=-4;
     if iLP==1
%     save lp
%     ' salvato ', keyboard
%     'mat_mix: da sistemare per profili termici non circolari ', keyboard
      Kmat_LP
     else
%      'mat_mix: da sistemare per profili termici non circolari ', keyboard
      Kmat_ve
     end
     toc
     
     ifp=ifp_salvo;

     igainshape=0;
     igainshapeT=0;

     %' dopo Kmat ', keyboard
     % lo salto perche` e` gia` fatto dentro a *_ci.m
     si2=[length(Pusa) length(Pusa)];
     if iLP==1
      KTem_p=KTempp(Pusa,Pusa,:);
      KTem_m=KTempm(Pusa,Pusa,:);
     else
      KTem_p=KTempp(Pusa,Pusa,:);
      KTem_m=KTempp(Pusa,Pusa,:);
      KTep_z=KTemzp(Pusa,Pusa,:);
      KTem_z=KTemzm(Pusa,Pusa,:);
     end

    else

       KTem_p=0;
       KTem_m=0;
       KTep_z=0;
       KTem_z=0;

    end

if iplan==0
KAp=KApsa;
KAm=KAmsa;
end
%'quia', keyboard

'% calcolo i coeff. di accop. della parte attiva'
%'quiattivo', keyboard

if ianti_gui==0
 yiN_ag=[];
 yiNd=[];
end
if ~exist('yiN')==1
 yiN=[];
 yiNd=[yiN yiN_ag];
end


if ( igau==0 | igau==5 | (igau==4 & length(yiNd)<=1) )

 if iplan==0
    if ~exist('KAp')
        disp(' errore 1 in mat_mix: manca KA area costante! ')
        pausak
    end
    if ~exist('KAp_ag')
     KAp_ag=KAp;
     KAm_ag=KAm;
    end
  end
elseif (abs(igau)==4 & length(yiNd)>1)

    disp(' calcolo i coeff. di accop. della parte attiva ')
%    keyboard
    %     sha=1;
    %     adis=xiN*kcav0;
    %     lxi=length(adis);
    igainshape=1;

    if ~exist('xiN')
        igainshape=2;
        global PGain PcGain PRef PcRef
        sN=abs(ismat(N));

%        if sN==1
        if max(sNt)==1
            adis=xiN*kcav;
            if length(find(shavet==4))>1
                ips=find(shavet==4 & real(nv(:,2))<2 & abs(imag(nv(:,2)))<.01 ...
                    & real(nv(:,2))>1);
                r_dr=radii.b(ips)/radii.a(ips); %attenzione: a b = y x
                sha=4;

                %         ' qui per set parametri ', keyboard
                iresK=0;
                bdis=adis*r_dr;
                pdi=pdi(1)*ones(size(bdis));
            end
%        elseif sN==2  %sN
        elseif max(sNt)==2  %sN


            sha=4;
            iresK=1;
            global PGain PcGain PRef PcRef
            roN=PcGain.xv;
            sP=size(PGain);
            Nr=sP(2);
            puP0=2:Nr;
            puP=1:Nr-1;
            xd=[]; yd=[]; pd=[]; vd=[];
            Ploca=PGain;
            Pm=(Ploca(:,puP)+Ploca(:,puP0))/2;
%            xd=[xd; Ploca(1,:)'];
%            yd=[yd; Ploca(2,:)'];
%            pd=[pd; Ploca(3,:)'];
            xd=[xd; Pm(1,:)'];
            yd=[yd; Pm(2,:)'];
            pd=[pd; Pm(3,:)'];
            val=Ploca(4,:)';
            dval=diff(val);
            dval=[dval; dval(end)];
            yiN=dval;


%            %           g0=PcGain.Glos;
%                      ' g0 ', keyboard

            Ploca=PRef;
            Pm=(Ploca(:,puP)+Ploca(:,puP0))/2;
%            xd=[xd; Ploca(1,:)'];
%            yd=[yd; Ploca(2,:)'];
%            pd=[pd; Ploca(3,:)'];
            xd=[xd; Pm(1,:)'];
            yd=[yd; Pm(2,:)'];
            pd=[pd; Pm(3,:)'];
%            xd=[xd; Ploca(1,:)'];
%            yd=[yd; Ploca(2,:)'];
%            pd=[pd; Ploca(3,:)'];
            val=Ploca(4,:)';
            dval=diff(val);
            dval=[dval; dval(end)];
            yiN_ag=dval*fatqw*NQW;

            xiT=xd+yd+pd;
            bdis=xiT*kcav;
            yiT=xd-yd+pd;
            adis=yiT*kcav;
            pdi=pd;


            %         ' spnmo in ,at_mix', keyboard
        end
        lxi=length(adis);

    else %xiN modo vecchio

        igainshape=1;

%        keyboard
        % parte vecchia
        fiat=find(iauto(:,1)==2);
        sha=shavet(fiat,1);
        shape=sha_fun(sha);
        rassi=radii.b(fiat,1)/radii.a(fiat,1);
        adis=xiN*kcav;
        bdis=xiN*kcav*rassi;
        pdi=ones(size(xiN))*radii.c(fiat,1);
        lxi=length(adis);

    end

     'Kmat attivo', keyboard
    ifp_salvo=ifp;

    if iLP==1
        %      'Kmat_LP', keyboard

        Kmat_LP
    else
        %      'Kmat_ve', keyboard
        Kmat_ve
    end
    igainshape=0;
    ifp=ifp_salvo;

    disp(' dopo calcolo i coeff. di accop. della parte attiva ')
%    keyboard

    if sha==1  %altri casi gia fatti in Kmat
        if iLP==1
           if sN>0
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
           end
           if sNr>0

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
                KAp_ag=KA_ag;
                KAm_ag=KA_ag;
               end
%'antigui', keyboard
            end

            KAp=KAp(Pusas,Pusas);
            KAm=KAp;
%                   'qui  LP', keyboard

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
            KAp=[KEEp segem*KEMp; segem*KMEp KMMp];
            KAm=[KEEp segem*KEMp; segem*KMEp KMMp];
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
                KAp_ag=[KEEp segem*KEMp; segem*KMEp KMMp];
                KAm_ag=[KEEp segem*KEMp; segem*KMEp KMMp];
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
            if ianti_gui==1
             KAp_ag=KAp_ag(Pusas,Pusas);
             KAm_ag=KAm_ag(Pusas,Pusas);
            end
            KAp=KAp(Pusas,Pusas);
            KAm=KAm(Pusas,Pusas);
            if iztm==1
                if ianti_gui==1
                 KAzp_ag=KAzp_ag(Pusas,Pusas);
                 KAzm_ag=KAzm_ag(Pusas,Pusas);
                end
                KAzp=KAzp(Pusas,Pusas);
                KAzm=KAzm(Pusas,Pusas);
            end
        else
            if ianti_gui==1
             KAp_ag=KAp_ag(Pusas,Pusas);
             KAm_ag=KAm_ag(Pusas,Pusas);
            end
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
 if ~exist('KAp')
  li=length(MKosp);
  for k=1:li
   du=MKosp{k};
   if prod(size(du))>0
    KAp=du(:,:,1);
    break
   end
  end 


 end
    map(KAp)
    title(' KAp ')
    pausak
   if iclo==1
%    'qui close ', keyboard
    hde=figure;
    close(hde-4:hde)
   end
end

% calcolo i coeff. di accop. dell'anisotropia e Idelta
%
%sha
%keyboard


%Kd=MKosp{6};
%' iredmat ', keyboard

if iredmat==1
    % Idelta=Idelta(Pusas,Pusas);
    % if iztm==1
    %  Ideltaz=Ideltaz(Pusas,Pusas);
    % end
    %' sav', keyboard
    if length(KAp)>length(Pusasav)
        KAp=KAp(Pusasav,Pusasav);
        KAm=KAm(Pusasav,Pusasav);
        if iztm==1
            KAzp=KAzp(Pusasav,Pusasav);
            KAzm=KAzm(Pusasav,Pusasav);
        end
    else
        ' in mat_mix ho commentato Pusasav !!!!!!!',
        if ifp==-10
            pausak
        end
    end

    sM=size(MKosp);
    MKospd=MKosp;
    MKosmd=MKosm;
    if iztm==1
        MKoszpd=MKoszp;
        MKoszmd=MKoszm;
    end
%    ' cosa succede?? ', keyboard
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
            for kj=1:kjM
                MKosp{ki}(Pus,Pus,kj)=MKospd{ki}(Pusasav,Pusasav,kj);
                MKosm{ki}(Pus,Pus,kj)=MKosmd{ki}(Pusasav,Pusasav,kj);
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
%' matmix qui ', keyboard
%    load Ks
%    Kosm{1}(:,:,2)=Ks1;
%    Kosp{1}(:,:,2)=Ks1;

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

si2=[length(Pusa) length(Pusa)];
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
%load sa

if ifp==-10
%' dopo mat mix', keyboard
end
%' dopo mat', keyboard
%' dopo mat', keyboard

%
%if length(inuo_bas)>0
%
% inuo_bas
%'nuo_bas', keyboard
%
% nuo_bas
%end

%if inume==0
 Kmat_any
%end
if icomp_anymat==1

%  'calcolo coef acc anisotropi', keyboard
  K_ANIP_lens
  if ifp==-10
   ' fine calcolo coef acc anisotropi', keyboard
  end 
end  

global imatm iret_fisso nosamix

if icomp_anymat==2
  if iret_fisso==1

   if length(imatm)==0
    kmat_vegra
    imatm=1;
    rad_anysav=rad_any;
      fipun=findstr(nomeFs,'.');
      radf=nomeFs(1:fipun-1);
      nosamix=['matmix-',radf];
     if exist([nosamix,'.mat'])==2
      raap=num2str(rand*1000);
      nosamix=[nosamix,'-',raap(1:3)];
     end
%     ' prima save   matmix', keyboard
    eval(['save ',nosamix,' Kplot1 Kplot2 KANp KANm KANzm KANzp rad_anysav'])
   else 
%    load matmix
%     ' prima   load matmix', keyboard

    eval(['load ',nosamix]);
   end
  else 
    kmat_vegra
  end 
%  'calcolo coef acc reticolo curvo', keyboard
%  kmat_vegra
  if ifp==-10
   ' fine calcolo coef acc ret. curvo', keyboard
  end 
end  

%global xprod yprod
%xprod=Tc;
%yprod=Tc;
% qui fine mat_mix', keyboard 