
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

  [vetd,aut]=eig_mio(N1,N2,nsol);
  vma=max(vetd);
  vet=vetd*diag(1./vma);

% DYG1=c./(aut*rr.*d*2);
%disp(' dopo di eig'), keyboard

%  [vetv,autv]=eig(N2,-N1);
% DYGv=c./(diag(autv)*rr.*d*2);

%disp(' dopo di eigs: verifica'), keyboard



 lam=1./aut;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PARTE 3: LEGAME CON LE CONDIZIONI DI SOGLIA %%%%%%%%%%%%%

  DYG=lam*c/(rr.*d*2);      %in 1/s in ampiezza

%  Gc0=DYG-DGp;
  Gc0=DYG;

  allam=(imag(Gc0)./real(Gc0));
%  Grlam=real(Gc0);
  Grlam=real(Gc0)-DGp;

  [Gv,ind]=sort(Grlam);
%  [Gv1,ind]=sort(Grlam);
  al=allam(ind);
%  To=Trasu*ones(1,length(ind));
%  V=vet(:,ind).*To;
  V=vet(:,ind);
%  V=diag(pes(Pus))*vet(:,ind);
%  'invpro', keyboard
  Lamv=lam(ind);
  Guama=1e17;

  igv=find(Gv>0 & Gv<Guama);
  if length(igv)<5
   igvd=find(Gv>0);
   if length(igvd)>5
    igv=igvd(1:5);
   else
    igv=igvd;
   end
%  elseif length(igv)>=31
%   fGmini=find(Gv(igv))<10*Gv(igv(1));
%   igv=igv(fGmini);
  end
if ifp==-11
 ' dopo autovalori: prima clear variabili'
 keyboard
end

  clear N1 N2 Mn11 Mn12 Mn21 Mn22 Mi11 Mi12 Mi21 Mi22

