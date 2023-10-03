%'acc0', keyboard
if exist('ikiaut')==1
 if length(ikiaut)==0
  clear ikiaut
 else
  if ikiaut==0
   clear ikiaut
  end
 end
end
if exist('ikiaut')==1
 if ikiaut==1
  isopro=1;
  izm=6;
  npuntik=2^(izm+1);
  kiniz=alimi;
  kfin=alim;
  [kiniz kfin]
  preca=.01;
  idis=0;
  iLPs=iLP;
  iLP=1;
  icalcola=1;
  acc01
  icalcola=1;
  iLP=iLPs;
  nk1=length(KK);
  isopro=0;

 if length(ifi0)==1
  icalcola=1;
  isopro=0;

 for iz=1:izm
  icalcola=1;
  ndk1=fix(npuntik/2^(iz+1));
  ndk=ndk1*2;

  if ifi0<=ndk1
   kiniz(iz)=KK(1);
   kfin(iz)=KK(ndk);
  elseif ifi0<=npuntik-ndk1
   kiniz(iz)=KK(ifi0-ndk1);
   kfin(iz)=KK(ifi0+ndk1);
  elseif ifi0>npuntik-ndk1
   kiniz(iz)=KK(ifi0-ndk1);
   kfin(iz)=alim;
  end
 end

 for iz=1:3
  clear Kp1
  ick=0;
  kiniz2=kiniz(iz);
  kfin1=kfin(iz);
  if iz==1
   kiniz1=alimi;
   kfin2=alim;
  else
   kiniz1=kiniz(iz-1);
   kfin2=kfin(iz-1);
  end
  Ks=KK;
  ick=0;
  for ik1=2:length(KK)
  k=abs(Ks(ik1));
   if (k<=kiniz2 & k>=kiniz1) | (k>=kfin1 & k<=kfin2)
    ick=ick+1;
    if ick<=2^(3-iz)-1
     KK(ik1)=-abs(KK(ik1));
    else
     ick=0;
    end
   else
    if ick~=0 & ik1<length(KK)
     while ick<2^(3-iz)-1
      ick=ick+1;
      KK(ik1)=-abs(KK(ik1));
     end
     ick=0;
    end
   end
  end
 end
  KK(1)=1e-8;
  fK=find(KK>=0);
  KK=KK(fK);

  dkf=(kfin(izm)-kiniz(izm))/4;
  kfin(izm+1)=kfin(izm)-dkf;
  kiniz(izm+1)=kiniz(izm)+dkf;

 for iz=4:izm
  clear Kp1
  ick=0;
  kinizi=kiniz(iz);
  kfini=kfin(iz);
  Ks=KK;
  for ik1=1:length(KK)-1
  k=Ks(ik1);
   if k>=kinizi & k<kfini
    ick=ick+1;
    Kp1(ick)=mean(Ks(ik1+[0 1]));
   end
  end
  Kt1=[Ks; Kp1'];
  KK=sort(Kt1);
 end
  if ifp>-2
   figure, subplot(211), plot(KK,'.-'), subplot(212), plot(diff(KK),'.-')
   pausak
  end
  'forem cotn',
  keyboard
  keyboard

  idis=2;
  ifp=ifpdum;
 end % length(ifi0)
  icalcola=1;
  acc01
  icalcola=1;
 end
else
 isopro=0;
 kiniz=alimi;
 kfin=alim;
 npuntik=nk1;
 acc01
end
