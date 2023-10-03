%icapol=input(' icapol (solo campo dominante = ');
%disp(' camv_new')
%keyboard





 if istrumix==0
  global fimaxi lfi fian0 Cug
  forma
  Cug.x=xcu;
%  'Cug', keyboard
  Cug.y=ycu;
  Cug.z=zcu;
 end



    a=a0ref*kcav0;
    avero=a0ref;
%    if length(xroJ)==0
%     disp(' vecchio modo no'), keyboard
%     Nx1=100;
%     Nx2=25;
%     if ifalso>=0
%      xx1=linspace(0,1*a,Nx1);
%      xx2=linspace(a+1e-2,amax,Nx2);
%     else
%      xx1=linspace(0,rax*a,Nx1);
%      xx2=linspace(rax+1e-2,rax*amax,Nx2);
%     end
%     xv=[xx1 xx2];
%     xvero=xv/a*avero;
%     x=xvero;
%
%    else
%     xvero=xroJ;
%     x=xvero;
%    end


    if length(xroJ)==1
     Nx=xroJ;
     plomax=2.5*max(max(axtot));
     plomay=2.5*max(max(aytot));
     ploma=max([plomax plomay]);
     xvero=linspace(0,ploma,Nx);
     x=xvero;
    else
     xvero=xroJ;
     x=xvero;
    end

%' xvero  ', keyboard

%    disp(' camful'), keyboard

    if length(lfi_inp)==0 | length(fimaxi)==0
     lfi_inp=51;
     fimaxi=2*pi;
    end

    Nx=length(x);
    npk=length(KK);


  lbv=length(mbv);
  fian=linspace(0,fimaxi,lfi_inp(1)+1);
%  fian=linspace(0,fimaxi,lfi_inp(1));
  if fimaxi<2*pi
   fian0=[fian(1:length(fian)-1) fian+pi];
  else
   fian0=fian;
  end

  fm0=cos(fian0);
  gm0=sin(fian0);

    xdx=[0 diff(xvero)].*xvero;
    defid=diff(fian0);
    defi=[defid defid(end)]';
    defi([1 end])=defi([1 end])/2;
    XP=xvero'*fm0;
    YP=xvero'*gm0;


    if idyn>=1
     if i2D==0
      fiGa=find(xvero>aiat);
      defi=2*pi;
     else
      borders
      Sat=ones(size(xvero'))*robor;
      XP1=xvero'*cos(fi);
      YP1=xvero'*sin(fi);
      Sbor=sqrt(XP1.^2+YP1.^2);
      fiGa=find(Sbor>Sat);
%      ' fiGA ', keyboard
     end
    end
    iLP=iLPr;
%    iLP=iLP1;

%' modangle', keyboard
%ipost=0
if ipost==0
    modangle
    mod_ff
    ipost=1;
end    

%global lp1 isavetutto lp2
%if length(isavetutto)==0
% isavetutto=0;
%end
%if isavetutto>0
%noms=['rad',num2str(isavetutto),'_',num2str(lp1),num2str(lp2)];
%clear Tstort Tstortp  MKoszmd MKoszpd MKoszd  MKospd MKosmd  Koszp Koszm  Kosp Kosm Kos Kosz Tstor Tstorb var
%clear Tstort Tstortp  MKoszmd MKoszpd MKoszd  MKospd MKosmd  Koszp Koszm  Kosp Kosm Kos Kosz Tstor Tstorb var
%clear DR ary ary_s Fi Tmirro Tstof Tstorbp Tdum Tcm Tc Pow Oo Odu Mns Mis
%eval(['save ',noms])
%end

    iLP=iLPr;


%' modangle dopo', keyboard

 al_shi=0;
% clear M2v Nazim Nrad tyE Lsov
 if length(if_only_out)==0
%  if ireset_int==-1
%   diss_new
%  else
%   diss_new
%   diss_nst
%   ' prima dissfun', keyboard
   eval(dissfun)
%   ' dopo dissfun', keyboard
%  end
 else
  if if_only_out==1
   diss_onl
  else
   diss_nst
  end
 end


