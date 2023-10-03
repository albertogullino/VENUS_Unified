if ~exist('Ofin')
 Ofin=1;
end 
%' calc_new per pu2 ', keyboard
%keyboard
if exist('Trcrit')==0
 Trcrit=eye(length(P)/2);
end

if exist('Trcritu')==1
 Trcrit=Trcritu{icpo};
end

ifpsave=ifp;
%ifp=-4;

if iraff==1
  freq=ze;
  fris=freq;
  [du,ifr]=min(abs(fre_camp-ze));
%  ifr=ifr+1;
%  ifr=1;
%'  ifast ', keyboard
  ifasasa=ifast;
  ifast=0;
  kcav=kcav0*(1+freq);
    %'scommenta !!!!'
    acc0_i
%    eval(['load ',nfsave])
 ifast=ifasasa;
  [du,ilam]=min(abs(imag(Lamv(Nd1vu))));
  iN=Nd1vu(ilam);
  fsol=iN;
%  Lami=Lamv(iN);
%  Ar=V(:,iN);

%    fim=find(Grlam>0 & Grlam<10*gg0);
%    Corre=abs(vet(:,fim)'*ADoutis);
%    if ifpsave>=-1
%     disp(' G alfa Corre ')
%     [Grlam(fim) allam(fim) Corre]
%    end
%    [du,fsold1]=min(abs(allam(fim)));
%    [du,fsold]=max(Corre);
%    if fsold~=fsold1
%     disp(' diversi')
%    end
%    fsol=fim(fsold);

    Ar=vet(:,fsol);   % autov all'interfaccia
%    Anout=Ar;
    glosout=2*Gv(fsol);
    [2*gg0 glosout]
   icsmax=sz(3);
   mmsav0=mmsav;   
%'qui campo ', keyboard
%'qui campo ', keyboard

 IOLD=0;   % e' gia' tutto fatto in inv_old
 if IOLD==1
      Apr=[Gas.*Ar; Ar];
      CCd=Apr;

      ics=1;
      Ab=CCd(pu2);
      Af=CCd(pu1);
      
      
      Afz(:,ics)=Ab.*Trs;
      Afzf(:,ics)=Af.*Trs;
%  '    passo', keyboard
%      Afz(:,ics)=Ab;
%      Afzf(:,ics)=Af;

      CCd=Tt*Apr;
      ics=ics+1;
      Ab=CCd(pu2);
      Af=CCd(pu1);
      Afz(:,ics)=Ab;
      Afzf(:,ics)=Af;

      Lam=Lamv(iN);
      Ttotal=(Mn+Lam*Mi);
      CCd=Ttotal*Apr;
      ics=ics+1;
      Ab=CCd(pu2);
      Af=CCd(pu1);
%      Afz(:,ics)=Ab.*Trd;
%      Afzf(:,ics)=Af.*Trd;
      Afz(:,ics)=Ab;
      Afzf(:,ics)=Af;


      sdu=size(Tmef);
      sTm=sdu(1:2);

%  sopra
      if icsfit>0
       for kf=1:icsfit
        ics=ics+1;
        CCd=reshape(Tmef(:,:,kf),sTm)*Apr;
        Ab=CCd(pu2);
        Af=CCd(pu1);
        Afz(:,ics)=Ab;
        Afzf(:,ics)=Af;
       end
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
        Afz(:,ics)=Ab;
        Afzf(:,ics)=Af;
       end
      end
      'fine OLD', keyboard
  else   %IOLD nuovo
   sAdu=size(Anz);
   Afz=reshape(Anz(:,1,:),sAdu([1 3]));
   Afzf=reshape(Anzf(:,1,:),sAdu([1 3]));
   
   ' Afz', keyboard
   Anout=(Afz(:,iFFsect)+Afzf(:,iFFsect));
   Anoute=Trcrit*Afz(:,1).*(Gas+1);
   Anout=Trcrit*Afz(:,1).*(Gas+1);
   Anouth=Trcrit*Afz(:,1).*(1-Gas);
%   'passo anout', keyboard
   AnQW=Afz(:,2)+Afzf(:,2);
   if iFFsect==1
    iFFa=3;
   else
    iFFa=1;
   end
   Anouta=Afz(:,iFFa)+Afzf(:,iFFa);   
  end  %IOLD
%disp(' fine autov '), keyboard
 clear Tb Tt Moi Mon Tmefb Tmef Mn Mi
 
   Ynor0=1./[Zve; Zvm];
   Ynor=Ynor0;
   
   for kn=1:numodiacc
    Ynor=[Ynor; Ynor0];
   end
   


   YD= diag(Ynor(Pusd));

%' Acoz raff', keyboard
   
   Acoz=[Afzf; Afz];
   if exist('Trcrit')==1
    Acoz(:,1)=[Trcrit*Afz(:,1).*Gas; Trcrit*Afz(:,1)];
   end
   Acoz_sav=Acoz;
   FF.Acoz=Acoz;
   duf=fmlstot(:,2);
   fidu=find(duf==0);
   duf(fidu)=1;
%   spesdu=Litot.*duf;
%   ' puend ', keyboard
   puout=find(iauto(:,1)==1);
%   puqw=1:find(iauto(:,1)==2)-1;
%   puend=1:find(iauto(:,1)==3)-2;
   puqw=puout:find(iauto(:,1)==2);
   puend=puout+1:find(iauto(:,1)==3)-1;
   duf0=duf;
   duf=zeros(length(dv),1);
   duf(puend)=duf0;
   spesdu=dv.*duf;
   lqw=sum(spesdu(puqw));
   lend=sum(spesdu(puend));

% critico
 if size(iauto,2)==3
   puouci=find(iauto(:,3)==-10);
   if length(puouci)==1
    puouc=1:puouci;
   else
    puouc=puouci(1):puouci(2);
   end
   lcri=sum(spesdu(puouc));
   FF.zi=[0 lcri lqw lend];   
   Nz(1)=nv(1,1);
   Nz(2)=nv(puouc(end),1);
   Nz(3)=rqw;
   Nz(4)=nv(end,1);
   FF.Nz=Nz;
  else 
   Nz(1)=nv(1,1);
   Nz(2)=rqw;
   Nz(3)=nv(end,1);
   FF.Nz=Nz;
   FF.zi=[0 lqw lend];   
  end 
%' end Acos', keyboard

else  %iraff

% if icampi==2
%  icsmax=ics;
% else
%  icsmax=2;
% end
%  icsmax=ics;

%' BAK ', keyboard, keyboard

% if length(if_only_out)==0
%   Afz=reshape(Azvet(:,pou(nso,fiAz),:,fiAz),sz([1 3]));
%   Afzf=reshape(Azvetf(:,pou(nso,fiAz),:,fiAz),sz([1 3]));
%   Anout=Afz(:,iFFsect)+Afzf(:,iFFsect);
%   AnQW=Afz(:,2)+Afzf(:,2);
%   ics=sz(3);
% else
%   Afz=reshape(Azvet(:,pou(nso,fiAz),fiAz),sz([1]),1);
%   Afzf=reshape(Azvetf(:,pou(nso,fiAz),fiAz),sz([1]),1);
%   Anout=Afz+Afzf;
%   AnQW=Afz+Afzf;
%   ics=1;
% end

%' verifica Amedio prima'
%keyboard
if length(if_only_out)==1
 if if_only_out==0
  if_only_out=[];
 end
end

 if length(if_only_out)==0
  if ifiez==0
%   ' calc_new 0 ', keyboard, keyboard
   ics=sz(3);
   icsmax=sz(3);
%   [Afz,sgu]=me_autn(Azvet,pou,nso,sz,Fint,ze,if_only_out,icsmax,ifpsave);
%   Afzf=me_autn(Azvetf,pou,nso,sz,Fint,ze,if_only_out,icsmax,ifpsave,sgu);
%   [Afz,sgu]=me_pro(Azvet,pou,nso,sz,Fint,ze,if_only_out,icsmax,ifpsave);
%save sutu Azvet pou nso sz Fint ze if_only_out icsmax ifpsave
   [Afz,sgu]=me_autu(Azvet,pou,nso,sz,Fint,ze,if_only_out,icsmax,ifpsave);
   Afzf=me_autu(Azvetf,pou,nso,sz,Fint,ze,if_only_out,icsmax,ifpsave,sgu);
%   Afzf=me_autu(Azvetf,pou,nso,sz,Fint,ze,if_only_out,icsmax,ifpsave);

%   Acoz(:,1)=[Gas.*Ar; Ar];
%   Acoz(:,2)=Aqw;
%   Acoz(:,3)=Aend;
%   ' calc_new 0 ', keyboard
   Ar=Afz(:,1)*0;
   Acoz=[Afzf; Afz]*0;
   if exist('Trcrit')==1
%    Acoz(:,1)=[Trcrit*Afzf(:,1); Trcrit*Afz(:,1)];
    Acoz(:,1)=[Trcrit*Afz(:,1).*Gas; Trcrit*Afz(:,1)];
%    'set Acoz', keyboard
  
%    Anout=Trcrit*Afz(:,1).*(Gas+1); 
   end
%   ' fine Acoz ', keyboard
   Acoz_sav=Acoz;
   Nz(1)=nv(1,1);
   Nz(2)=rqw;
   Nz(3)=nv(end,1);
   FF.Nz=Nz;
   FF.Acoz=Acoz;
   duf=fmlstot(:,2);
   fidu=find(duf==0);
   duf(fidu)=1;
%   spesdu=Litot.*duf;
%   ' puend ', keyboard
   puout=find(iauto(:,1)==1);
%   puqw=1:find(iauto(:,1)==2)-1;
%   puend=1:find(iauto(:,1)==3)-2;
   puqw=puout:find(iauto(:,1)==2);
   puend=puout+1:find(iauto(:,1)==3)-1;
   duf0=duf;
   duf=zeros(length(dv),1);
   duf(puend)=duf0;
   spesdu=dv.*duf;
   lqw=sum(spesdu(puqw));
   lend=sum(spesdu(puend));

   FF.zi=[0 lqw lend];
%   ' calc_new 0 ', keyboard
  elseif ifiez==1
%   ' calc_new 1 ', keyboard
   Afz(:,1)=Gas.*Ar;
   Afzf(:,1)=Ar;
   Afz(:,2)=Aqw(pu2);
   Afzf(:,2)=Aqw(pu1);
   Afz(:,2)=Aend(pu2);
   Afzf(:,2)=Aend(pu1);
  elseif ifiez==2
%   ' calc_new 2 ', keyboard
   [du,ilam]=min(abs(imag(Lamv(Nd1vu))));
   iN=Nd1vu(ilam);
   Lami=Lamv(iN);
   ggg1=real(Lami)*c/(rr.*d*2);
   Ar=V(:,iN);
   A0=[Gas.*Ar; Ar];
   Aqw=Tt*A0;
   Aqw1=(Mons+Mois*Lami)*Aqw;
   Aend=Tbfield*Aqw1;

   Afz(:,1)=Gas.*Ar;
   Afzf(:,1)=Ar;
   Afz(:,2)=Aqw(pu2);
   Afzf(:,2)=Aqw(pu1);
   Afz(:,2)=Aend(pu2);
   Afzf(:,2)=Aend(pu1);
  end


%   Anout=Trcrit*(Afz(:,iFFsect)+Afzf(:,iFFsect));

   Anout=(Afz(:,iFFsect)+Afzf(:,iFFsect));
   Anoute=Trcrit*Afz(:,1).*(Gas+1).*Ofin;
   Anout=Trcrit*Afz(:,1).*(Gas+1).*Ofin;
     iprima=0;
%     iprima=1;
%     'ipriomoa', keyboard
     if iprima==1
      Anout=[Afzf(:,iFFsect); Afz(:,iFFsect)];
      Anouts=[Afzf(:,iFFsect); Afz(:,iFFsect)];
%      Afzf_s=Afzf;
%      Afz_s=Afzc;
      load gaussOo4
  %    Oo=inv(OoScala);
      Oo=(OoScala);
      Anouted=Oo*Anout;
      Afzf(:,iFFsect)=Anouted(1:end/2);
      Afz(:,iFFsect)=Anouted(end/2+1:end);
      Anoute=Anouted(1:end/2)+Anouted(end/2+1:end);
      Anout=Anoute;
     end 

   Ynore0=1./[Zve];
   Ynorm0=1./[Zvm];
   Ynore=Ynore0;
   Ynorm=Ynorm0;
   
   for kn=1:numodiacc
    Ynore=[Ynore; Ynore0];
    Ynorm=[Ynorm; Ynorm0];
   end
   Ynor=[Ynore; Ynorm];

%   'Ybor', keyboard
   if ipolar==0
    Ynor=[Ynor Ynor];
   end

   YD= diag(Ynor(Pusd));


   Anouth= YD*Trcrit*Afz(:,1).*(1-Gas).*Ofin;
%   'passo anout', keyboard
   AnQW=Afz(:,2)+Afzf(:,2);
   if iFFsect==1
    iFFa=3;
   else
    iFFa=1;
   end
   Anouta=Afz(:,iFFa)+Afzf(:,iFFa);

 else
   ics=1;
   icsmax=ics;
   [Afz,sgu]=me_autn(Azvet,pou,nso,sz,Fint,ze,if_only_out,icsmax,ifpsave);
  % [Afz,sgu]=me_pro(Azvet,pou,nso,sz,Fint,ze,if_only_out,icsmax,ifpsave);
   Afzf=me_autn(Azvetf,pou,nso,sz,Fint,ze,if_only_out,icsmax,ifpsave,sgu);

   Anout=Trcrit*(Afz+Afzf);
   Anouta=Afz+Afzf;
   AnQW=Afz+Afzf;

 end
 
% ' Anotut in calcnew', keyboard
    mmsav0=mm;
    mmsav0=mmsav;
    if mm~=mmsav0 & iLP==0
%    if mm~=mmsav
     if itetm==1
      Anout=[Anout; zeros(size(Anout))];
      AnQW=[AnQW; zeros(size(AnQW))];
     else
      Anout=[zeros(size(Anout)); Anout];
      AnQW=[zeros(size(AnQW)); AnQW];
     end
   end

%' verifica Amedio dopo'
%keyboard

 Anff=Anout;
 AnQWs=AnQW;


   if iant==0
    Anout=Anoute.*pes(Pus);
    AnQW=AnQW.*pes(Pus);
   end
%   'A3', keyboard


end

faEzsav=fatEz;

ICSM=icsmax;

if imazver>1
 ICSM=1;
end

 for kmo=1:ICSM
%   Andu=Afz(:,kmo);
  if kmo==1
   Andu=Anoute;
   AnduH=Anouth;
  else
   Andu=(Afz(:,kmo)+Afzf(:,kmo));
   AnduH=YD*(Afz(:,kmo)-Afzf(:,kmo));
  end
  if length(Andu)<3
   Andu=Anoute;
   AnduH=Anouth;
  end
%   ' calc_new', keyboard

    if mm~=mmsav0 & iLP==0
%    if mm~=mmsav
     if itetm==1
      Andu=[Andu; zeros(size(Andu))];
      AnduH=[AnduH; zeros(size(AnduH))];
     else
      Andu=[zeros(size(Andu)); Andu];
      AnduH=[zeros(size(AnduH)); AnduH];
     end
    end
       if ired_ret==1
        A1=zeros(ldap(end),1);
        A1(Pured1)=Andu;
        A2=zeros(ldap(end),1);
        A2(Pured1)=AnduH;
        Andu=A1;
        AnduH=A2;
       end

%  plot FF
   if kmo==1
    fatEz=faEzsav;
   else
    fatEz=1;
   end

%kmo
%   'cam_val prima',   keyboard

   cam_val_arm;
   %'cam_val_arm dopo',   keyboard


     Pol.Ex=0;
     Pol.Ey=0;

   if iskim==1, break, end
%   ' sono qui calc_new', keyboard
%iallef=iall;
iallef=1;
%   if kmo==iFFsect & iFF==10
   if kmo==iFFsect & iFF==1
%     AnFF=Andu;
%     plo_ff
  if iallef==1
     p_ff_new
  else
   ' no FF', keyboard
    Ef=0;
    X=0;
    Y=0;
    Efx=0;
    Efy=0;
  end  
   end


   if kmo==iFFsect
%     Anduf=Trcrit*Afzf(:,kmo);
     
%     if mm~=mmsav0
%      if itetm==1
%       Anduf=[Anduf; zeros(size(Anduf))];
%      else
%       Anduf=[zeros(size(Anduf)); Anduf];
%      end
%     end
%        if ired_ret==1
%         A1=zeros(ldap(end),1);
%         A1(Pured1)=Anduf;
%         Anduf=A1;
%        end

%     AnFF=Andu+Anduf;
%     AnFFH=Andu-Anduf;
     AnFF=Andu;
     AnFFH=AnduH;
     
    E2x=E2xd;
    E2y=E2yd;
    if isalva==1
     Pol.Ex=Exdu;
     Pol.Ey=Eydu;
     isalva=0;
    end
%    mod_cha
%    disp(' sezione uscita '), keyboard
   end

   
   if i2D==3
    Etx(:,:,kmo)=E2xd;
    Ety(:,:,kmo)=E2yd;
%    Etx(:,:,kmo)=Exdu;
%    Ety(:,:,kmo)=Eydu;    
    Htx(:,:,kmo)=Hxdu;
    Hty(:,:,kmo)=Hydu;        
%    'scrivo Ht', pausak
    if iztm==1
     Etz(:,:,kmo)=E2zd;
    else
     Etz(:,:,kmo)=E2yd*0;
    end
    
   else
    Etb(:,kmo)=Etot;
    Etx(:,kmo)=E2xd;
    Ety(:,kmo)=E2yd;
    if iztm==1
     Etz(:,kmo)=E2zd;
    else
     Etz(:,kmo)=E2yd*0;
    end
   end
   if kmo==2
      Andu=Afz(:,kmo)+Afzf(:,kmo);
      if mm~=mmsav0 & iLP==0
%      if mm~=mmsav
       if itetm==1
        Andu=[Andu; zeros(size(Andu))];
       else
        Andu=[zeros(size(Andu)); Andu];
       end
     end
    if ired_ret==1
     A1=zeros(ldap(end),1);
     A1(Pured1)=Andu;
     Andu=A1;
    end
      cam_val;
      Etcsum=Etot;
   end

   if idyn>=1
    Poz=(Exdu.*conj(Hydu)-Eydu.*conj(Hxdu));
    Poi(:,:,kmo)=Poz;    
    %' Poi ver', keyboard
   end
%' ora salto questa parte', keyboard
   if idyn>=1
%   if idyn>=100
%    if kmo==2
%     Andu=Afz(:,kmo)+Afzf(:,kmo);
%    else
     Andu=Afz(:,kmo);
        if ired_ret==1
         A1=zeros(ldap(end),1);
         A1(Pured1)=Andu;
         Andu=A1;
        end
%    end
    cam_val;
    Etotb=Etot;
    Andu=Afzf(:,kmo);
        if ired_ret==1
         A1=zeros(ldap(end),1);
         A1(Pured1)=Andu;
         Andu=A1;
        end
    cam_val;
    if i2D==3
     Etf(:,:,kmo)=Etot;
     Etb(:,:,kmo)=Etotb;
    else
     Etf(:,kmo)=Etot;
     Etb(:,kmo)=Etotb;
    end
   end
 end

ifp=ifpsave;
fatEz=faEzsav;

