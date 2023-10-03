  iaut_new=1;
  %iaut_new=0;
  if ired_ret==1
    iaut_new=2;
  end
%  igam=1
  igam=0;
%' inv_old ', keyboard
%' inv_old ',
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
   Mns=Mn;
   Mis=Mi; 

  if (igam==1 | ifr==1) & pol==pvet(1)
%   ' in inv_old', keyboard
   Gadu=Gad;
   dRi=diag(Gadu)*Mn12-Mn22;
   dPi=Mn21-diag(Gadu)*Mn11;
   Geqm=inv(dRi)*dPi;  % riflessione di tutta la struttura da sopra
   
%   Gadu=Ga2;
%   Ap=Mn21*Gadu+Mn22;
%   Bp=Mn11*Gadu+Mn12;
%   Geqm=Bp*INv(Ap);
   

     Geq=abs(diag(Geqm));

%    if  ifp==-10 %& ifiez>0
   if max(Geq)>.1
     if ifp==-10
      hpro=figure; plot(Geq), 
      pausak
     end   
%   if max(Geq)>1 | Ps.igacrit>=1
%     if ifp==-10
%      hpro=figure; plot(Geq), 
%     keyboard
%     end
     %pausak,      keyboard
%     close
   end 
%  King=KKt(Pusas(1:ldapuu(2)));
  if length(Pusas)>130
   istop=0;
  else
   istop=1;
  end

   if igam==1
    keyboard
   end
  end

%Ga1=diag(Gas);
%Ga2=diag(Gad);

  N1=Mn22-Ga2*Mn11*Ga1+Mn21*Ga1-Ga2*Mn12;
  N2=Mi22-Ga2*Mi11*Ga1+Mi21*Ga1-Ga2*Mi12;

%' ENNE', keyboard
%  N2a=Mi22;
%  N2b=-Ga2*Mi11*Ga1;
%  N2c=Mi21*Ga1;de
%  N2d=-Ga2*Mi12;
%  N2=Mi22-Ga2*Mi11*Ga1+Mi21*Ga1-Ga2*Mi12;
%  N2s=N2a+N2b+N2c+N2d;

%   N1i=-inv(N1);
%   Mtot=N1i*N2;
%  [A1,B1]=eig(N2,-N1);
%' inv_old  6666', keyboard
%' inv_old ', keyboard
N1s=N1;
N2s=N2;
 if ifp==-10 & ifr==1 & igam==2
  figure, plot(abs(diag(N1s))), pausak
 end

%' inv_old ', keyboard

  nsol=round(length(N1)/10);
  if nsol>10
   nsol=10;
  end
  nsol=0;  %tutte
%disp(' prima di eig'), keyboard

global PUp1 PUp2 polsub
  polsub=pol;
  if ifr==1 & pol==-1
   clear PUp1
   clear global PUp1
  end
  if ifr==1 & pol==1
   clear PUp2
   clear global PUp2
  end

  if iaut_new==1
   [vetd,aut]=eig_mio(N1,N2,nsol);
%   'new'
  elseif iaut_new==2
   [vetd,autd]=eig(N2,-N1);
   aut=diag(autd);
  else
%   N1i=-inv(N1);
%   Mtot=N1i*N2;
	%    [A1,B1]=eig(Mtot);
%    [A2,B2]=eig(Mtot,'nobalance');
    N1i=-inv(N1);
    Mtot=N1i*N2;
    [vetd,autd]=eig(Mtot,'nobalance');
    aut=diag(autd);
%    N1i=-inv(N1s);
%    Mtot=N1i*N2s;
%    [vetd2,autd2]=eig(Mtot,'nobalance');
%    aut2=diag(autd2);    
  end
%' inv_old ', keyboard

  vma=max(vetd);
  vet=vetd*diag(abs(vma)./vma);
  vet=vetd;

% DYG1=c./(aut*rr.*d*2);
%disp(' dopo di eig'), keyboard

%  [vetv,autv]=eig(N2,-N1);
% DYGv=c./(diag(autv)*rr.*d*2);

%if ifp==-10
%disp(' dopo di eigs: verifica'), pausak
%end


 lam=1./aut;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PARTE 3: LEGAME CON LE CONDIZIONI DI SOGLIA %%%%%%%%%%%%%

%  DYG=lam*c/(rr.*d*2);      %in 1/s in ampiezza
%  DYGo=lam*c/(rr.*d*2);      %in 1/s in ampiezza

%  DYG=lam*c*rr/(ra^2.*d*2);      %in 1/s in ampiezza
%  DYG=lam*rr/(real(ra).*d*2*1e2);      %in 1/s in ampiezza
  DYG=lam*rr/(real(ra).*d*2*1e2);      %in 1/s in ampiezza

%  DYG=lam*c/(rr.*d*2);      %in 1/s in ampiezza

%lamf=P*n/(2*na*La*1e-4);
  Gc0=DYG+DGp;
%  Gc0=DYG;
%'Cg0 cont', keyboard

  allam=(imag(Gc0)./real(Gc0));
  Grlam=real(Gc0);
%  Grlam=real(Gc0)-DGp;


  [Gv,ind]=sort(Grlam);
%  [Gv1,ind]=sort(Grlam);
%'Dyg', keyboard
%'Dyg', keyboard

  al=allam(ind);
  
%  fiv=find(Gv>0);
%  Gva=Gv(fiv);
%  ava=al(fiv);

%  To=Trasu*ones(1,length(ind));
%  V=vet(:,ind).*To;
  V=vet(:,ind);
%  V=diag(pes(Pus))*vet(:,ind);
%  'invpro', keyboard
  Lamv=lam(ind);
  Guama=2e15;
%  Guama=1e17;

  igv=find(Gv>0 & Gv<Guama);
  if length(igv)<5
   igvd=find(Gv>0);
   if length(igvd)>5
    igv=igvd(1:5);
   else
    igv=igvd;
   end
  else
%  elseif length(igv)>=fix(1.5*nmasce)
   fGmini=find(Gv(igv)<1000*mean(Gv(igv(1:3))));
	   igv=igv(fGmini);
  end
% ' dopo autovalori: prima clear variabili'
% keyboard
if ifp==-111
 ' dopo autovalori: prima clear variabili'
 keyboard
end

%if ifr==21
%' fine invmix', keyboard   
%end

%' BAK ', keyboard, keyboard

  clear N1 N2 Mn11 Mn12 Mn21 Mn22 Mi11 Mi12 Mi21 Mi22


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
     keyboard
     end
if ifr==8 | ifr==1
% 'controllo merda', keyboard
end
if ifp==-11 | ifp==0
 disp(' dopo inv_mix ')
 keyboard
end
   if isopro>1
    return
   end

  if iraffi==1
   %'iraffi', keyboard
    [du,ilam]=min(abs(imag(Lamv(Nd1vu))));
    iN=Nd1vu(ilam);
    Nd1vu=iN;
  end  

    iNc=0;

    for iN=Nd1vu'
      iNc=iNc+1;
      Ar=V(:,iN);
      Lam=Lamv(iN);
      
      Apr=[Ga1*Ar; Ar];
      CCd=Apr;
      ics=1;
      Af=CCd(pu1);
      Ab=CCd(pu2);
       Afm=Af;
       Abm=Ab;      
     if icriso==100           %non devo farlo, perche' e' gia' in calc_new!!!!

       Ab=Trcrit*Ar;     
       Af=Ab.*Gas;     
     end
     %'gacrit verifica', keyboard
      Anz(:,iNc,ics)=Ab;
      Anzf(:,iNc,ics)=Af;

%      Ts{ifr}=Tt;

      CCd=Tt*Apr;
      ics=ics+1;
      Ab=CCd(pu2);
      Af=CCd(pu1);
      Anz(:,iNc,ics)=Ab;
      Anzf(:,iNc,ics)=Af;
      Abqw=Ab;

      Ttotal=(Mn+Lam*Mi);
      CCd=Ttotal*Apr;
      ics=ics+1;
      Ab=CCd(pu2);
      Af=CCd(pu1);
%      Anz(:,iNc,ics)=Ab.*Trd;
%      Anzf(:,iNc,ics)=Af.*Trd;
      Anz(:,iNc,ics)=Ab;
      Anzf(:,iNc,ics)=Af;




      if exist('Tmef')
       sdu=size(Tmef);
       sTm=sdu(1:2);
      end


%  sopra
%disp('icsfit'), keyboard

      if icsfit>0 & icriso==0
%       for kf=1:icsfit
%        ics=ics+1;
%        CCd=reshape(Tmef(:,:,kf),sTm)*Apr;
%        Ab=CCd(pu2);
%        Af=CCd(pu1);
%        Anz(:,iNc,ics)=Ab;
%        Anzf(:,iNc,ics)=Af;
%       end
        Tm=1;
         for kf=1:icsfit
          Tm=Tmef(:,:,kf)*Tm;
         end
        CCd=Tm*Apr; 
        Abm=CCd(pu2);
        Afm=CCd(pu1);
      end

      if icsfit>0
        ics=ics+1;
        Anz(:,iNc,ics)=Abm;
        Anzf(:,iNc,ics)=Afm;
        istmet=ics;
      end	
      if exist('Tmefb')
       sdu=size(Tmefb);
       sTm=sdu(1:2);
      end
%      ' da fare', keyboard
    if idyn>0
       TQW=(Mon+Lam*Moi);
       AprQW=TQW*Tt*Apr;
        ics=ics+1;
        CCd=AprQW;
        Ab=CCd(pu2);
        Af=CCd(pu1);
        Anz(:,iNc,ics)=Ab;
        Anzf(:,iNc,ics)=Af;      
     end   
%  sotto
      if icsfib>0
      ' da fare', keyboard
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

    if exist('isaMo')
     if isaMo==1
      Mons=Mon;
      Mois=Moi;
     end
    end
    if ifiez~=2 & iraffi==0
%    if ifiez~=2
     clear Tb Tt
    else
     Mons=Mon;
     Mois=Moi;
    end
%    clear Moi Mon Tmefb Tmef Mn Mi

%'inv', keyboard