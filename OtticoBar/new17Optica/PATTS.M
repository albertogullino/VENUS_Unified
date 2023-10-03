clear
close all
load lav
ifp=-10;

    if ( igau==0 | (igau==4 & length(yiN)<=1) )

     if ~exist('KAp')
      disp(' errore 1 in mat_mix: manca KA area costante! ')
      pausak
     end

    elseif (abs(igau)==4 & length(yiN)>1)

disp(' calcolo i coeff. di accop. della parte attiva ')
keyboard
%     sha=1;
%     adis=xiN*kcav0;
%     lxi=length(adis);
     igainshape=1;

     if ~exist('xiN')
      igainshape=2;
       global PGain PcGain PRef PcRef
       sN=size(N);

       if length(sN)==1
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
       elseif length(sN)==2  %sN


          sha=4;
          iresK=1;
          global PGain PcGain PRef PcRef
          roN=PcGain.xv;
          sP=size(Pterm{1});
          Nr=sP(2);
          xd=[]; yd=[]; pd=[]; vd=[];
          Ploca=PGain;
           xd=[xd; Ploca(1,:)'];
           yd=[yd; Ploca(2,:)'];
           pd=[pd; Ploca(3,:)'];
           val=Ploca(4,:)';
           dval=diff(val);
           dval=[dval; dval(end)];
           yiN=dval;
%           g0=PcGain.Glos;
%           ' g0 ', keyboard

          Ploca=PRef;
           xd=[xd; Ploca(1,:)'];
           yd=[yd; Ploca(2,:)'];
           pd=[pd; Ploca(3,:)'];
           val=Ploca(4,:)';
           dval=diff(val);
           dval=[dval; dval(end)];
           yiN_ag=dval;

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

      keyboard
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
    end

% 'Kmat attivo', keyboard
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
keyboard
