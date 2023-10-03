  s=size(Mn);
  l1=s(1)/2;
  l2=s(1)/2+1;
  l3=s(1);
  Mn11=Mn(1:l1,1:l1);
  Mn12=Mn(1:l1,l2:l3);
  Mn21=Mn(l2:l3,1:l1);
  Mn22=Mn(l2:l3,l2:l3);
  Mi11=Mi(1:l1,1:l1);
  Mi12=Mi(1:l1,l2:l3);
  Mi21=Mi(l2:l3,1:l1);
  Mi22=Mi(l2:l3,l2:l3);


  N1=Mn22-Ga2*Mn11*Ga1+Mn21*Ga1-Ga2*Mn12;
  N2=Mi22-Ga2*Mi11*Ga1+Mi21*Ga1-Ga2*Mi12;


  nsol=round(length(N1)/10);
  if nsol>10
   nsol=10;
  end
  nsol=0;  %tutte
%disp(' prima di eig'), keyboard

 if iaut_d==0
  [vet,aut]=eig_mio(N1,N2,nsol);
 else
  [vet,aut]=eig_mio(N1.',N2.',nsol);
 end



  lam=1./aut;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PARTE 3: LEGAME CON LE CONDIZIONI DI SOGLIA %%%%%%%%%%%%%

  DYG=lam*c/(rr.*d*2);      %in 1/s in ampiezza

%  Gc0=DYG-DGp;
  Gc0=DYG;

  allam=(imag(Gc0)./real(Gc0));
%  Grlam=real(Gc0);
  Grlam=real(Gc0)-DGp;

  [Gv,ind]=sort(Grlam);
  al=allam(ind);
  V=vet(:,ind);
  Lamv=lam(ind);
  Guama=1e17;

  igv=find(Gv>0 & Gv<Guama);

%  Lam_rig=Lamv(igv);

  if length(igv)<5
   igvd=find(Gv>0);
   if length(igvd)>5
    igv=igvd(1:5);
   else
    igv=igvd;
   end
  end
if ifp==-11
 ' dopo autovalori: prima clear variabili'
 keyboard
end

%  clear N1 N2 Mn11 Mn12 Mn21 Mn22 Mi11 Mi12 Mi21 Mi22




   if isopro==1
    return
   end

     [alphavv,iav]=sort(al(igv));
     nmodv=[1:length(iav)];
     Nd1vu=igv(iav);
     Gvav=Gv(Nd1vu);


     if ifp>=-1
     disp(' fr , Gvav, alphav')
     [ Gvav alphavv]
     end

if ifp==-11 | ifp==0
 disp(' dopo inv_mix ')
 keyboard
end
   if isopro>1
    return
   end

  if iraffi==0

    iNc=0;

    for iN=Nd1vu'
      iNc=iNc+1;
      Ar=V(:,iN);
%      ' Lamris '
      Lam=Lamv(iN);
      Lam_rig(iNc)=Lam;
%      pausak
      Apr=[Gas.*Ar; Ar];
      CCd=Apr;
      ics=1;
      Ab=CCd(pu2);
      Af=CCd(pu1);
%      Anz(:,iNc,ics)=Ab.*Trs;
%      Anzf(:,iNc,ics)=Af.*Trs;
      Anz(:,iNc,ics)=Ab;
      Anzf(:,iNc,ics)=Af;

      CCd=Tt*Apr;
      ics=ics+1;
      Ab=CCd(pu2);
      Af=CCd(pu1);
      Anz(:,iNc,ics)=Ab;
      Anzf(:,iNc,ics)=Af;

      Ttotal=(Mn+Lam*Mi);
      CCd=Ttotal*Apr;
      ics=ics+1;
      Ab=CCd(pu2);
      Af=CCd(pu1);
%      Anz(:,iNc,ics)=Ab.*Trd;
%      Anzf(:,iNc,ics)=Af.*Trd;
      Anz(:,iNc,ics)=Ab;
      Anzf(:,iNc,ics)=Af;


%' BAK ', keyboard, keyboard

      if exist('Tmef')
       sdu=size(Tmef);
       sTm=sdu(1:2);
      end


%  sopra
%disp('icsfit'), keyboard

      if icsfit>0
       for kf=1:icsfit
        ics=ics+1;
        CCd=reshape(Tmef(:,:,kf),sTm)*Apr;
        Ab=CCd(pu2);
        Af=CCd(pu1);
        Anz(:,iNc,ics)=Ab;
        Anzf(:,iNc,ics)=Af;
       end
      end

      if exist('Tmefb')
       sdu=size(Tmefb);
       sTm=sdu(1:2);
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
        Anz(:,iNc,ics)=Ab;
        Anzf(:,iNc,ics)=Af;
       end
      end

    end   %iNc

%    clear Tb Tt Moi Mon Tmefb Tmef Mn Mi

   end % iraffi

