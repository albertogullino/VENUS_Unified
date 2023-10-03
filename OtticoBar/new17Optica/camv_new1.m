%icapol=input(' icapol (solo campo dominante = ');
%disp(' camv_new')
%keyboard
if ~exist('ipost')==1
 ipost=0;
end


if ipost==1
 dissfun='diss_nst2';
% dissfun='diss_nst1';
% dissfun='diss_nst0';
%  dissfun='diss_nst7';
%  dissfun='diss_nst5';
 dissfun='diss_nst4_wol';
else
%dissfun='diss_nst';
 dissfun='diss_nst1';
 dissfun='diss_nst1_wol';
 dissfun='diss_nst2_wol';
 dissfun='diss_nst3_wol';
 dissfun='diss_nst4_wol';
% dissfun='diss_nst3';
% dissfun='diss_nst2';
% dissfun='diss_nst7';
% dissfun='diss_nst0';
 global irifo
  if irifo==1
   dissfun='diss_nst1';
  end

end
%dissfun='diss_new';
%dissfun='diss_fie';
%dissfun='diss_las';
%dissfun='diss_nst1';

icapol=1;
if length(icapol)==0
 icapol=0;
end
if idyn>=1
 icampi=2;
end
if iLP==1
 icapol=0;
end


%if ifps<=-2
 if fix(nube/2)-nube/2==0
  desi=0;
 else
  desi=1;
 end
%end

if icampi>0
  if  length(desi)==0
   desi=input(' plot a sinistra [1] o destra [0] = ');
   if length(desi)==0
    desi=0;
   end
  end
  enlar=0;
  nsub2=1;
  enlar=300;
  nsub2=2;
    ved=version;
    ve=ved(1);
  if iLP==0
   if desi==1
    sinc=+1;
    pograp=[ 40 70 450 600];
    pogram=[ 40 30 450 600];

    if str2num(ve(1))==7
     pograp=[33   200   755   772];
     pogram=[33   104   755   772];
     pograp=[20   200   500   550];
     pogram=[20   100   500   550];
    end 
    else
    sinc=-1;
    pograp=[ 550 70 450 600];
    pogram=[ 550 30 450 600];

    if str2num(ve(1))==7
     pograp=[550   200   755   772];
     pogram=[550   104   755   772];
    end 
   end
  else
   if desi==1
    sinc=+1;
    pograp=[0 150 450 500];
    pogram=[0 50 450 500];
    else
    sinc=-1;
    pograp=[550 150 450 500];
    pogram=[550 50 450 500];
   end
  end

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

%' xvero  '
%keyboard

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


end   %icampi
%' modangle dopo', keyboard

 al_shi=0;
% clear M2v Nazim Nrad tyE Lsov
 if length(if_only_out)==0
%  if ireset_int==-1
%   diss_new
%  else
%   diss_new
%   diss_nst
   eval(dissfun)
%  end
 else
  if if_only_out==1
   diss_onl
  else
   diss_nst
  end
 end


