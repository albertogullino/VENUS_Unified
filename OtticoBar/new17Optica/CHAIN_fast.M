 if ifr==1 & pol==pvet(1)
%  ' entro in chain_fast',   keyboard
 end 
%ifp=-11
fmlsin=[];
Tstor=[];
%ifp=-11
%ifpss=ifp;
%ifp=-4;
if ifp==-12
' entro in chain_i', keyboard
end
%ssh=size(shai);
%if ssh(2)>1
% shailoc=reshape(shai,1,prod(ssh));
%else
% shailoc=shai;
%end
clear Ide Tmi
if ifiez==1
' NB: field along z only for Temp=0'
end

shailoc=shai(1,:);

if inuo>=1
% fisem=find(ai-bi~=0);
 sni=sum(ni);
 fisem=find(ai(1,:)~=0 & Li'<1);
% fisem=find(ai(1,:)~=0);
% asem=zeros(1,length(ai));
 asem=zeros(1,length(ai(1,:)));
 asem(fisem)=ai(1,fisem);
% ilaymem=find(diff([asem 0])<0 & shailoc>1);
% ilaymem=find(diff([asem 0])<0 );
% ilave=find(asem~=0 & shailoc>1);
 ilave=find(asem~=0 & shailoc>0 | fmlsi(:,2)'>0);
 ilaymem=zeros(1,length(ai));
 ilaymem(ilave)=1;
%' modificata ilaymem !!!!!!!!!!!!'
 ilaymem(1:end)=1;
 if inuo==1
  if pol==pvet(1)
   layskip=zeros(1,length(ai));
   ilaski=0;
  else
   filay=find(asem~=0 & shailoc>1);
   layskip=ones(1,length(ai));
   layskip(filay)=0;
   ilaski=1;
  end
 else
   layskip=zeros(1,length(ai));
   ilaski=0;
 end

else

 ilaski=0;
% layskip=zeros(size(ai));
 layskip=zeros(1,length(ai));
 ilaymem=zeros(1,length(ai));
end

 ilayfastd=zeros(1,length(ai));
 if ifastsa==0
  fisem=find(ai(1,:)~=0);
%  asem=zeros(1,length(ai));
  asem=zeros(1,length(ai(1,:)));
  asem(fisem)=ai(1,fisem);
  ilave=find(asem~=0 & shailoc>0 & Li'<1);
%  ilave=find(asem~=0 & shailoc>0);
%  ilave=find(asem~=0 & shailoc>0 & abs(imag(sni))<.1);
  ilayfastd(ilave)=1;
  ilayfast=ilayfastd;
 elseif ifastsa==1
  fisem=find(ai(1,:)~=0);
%  asem=zeros(1,length(ai));
  asem=zeros(1,length(ai(1,:)));
  asem(fisem)=ai(1,fisem);
  ilave=find(asem~=0 & shailoc>0);
  ilayfast=ilayfastd;
  ilayfastd(ilave)=1;
  if itop==1
   ilayfastto=ilayfastd;
  else
   ilayfastbo=ilayfastd;
  end
 else
  ilayfast=ilayfastd;
 end
%' ilay mem', keyboard
%sMAT=length(find(ilaymem==1))*length(Pust)^2*16;
%Mb=2^20;
%
%if iTsav==1 & ifr==1 & pol==pvet(1) & icredir==1 & sMAT>Mb*300
%  ck=clock;
%  Dire=num2str(floor(prod(ck( find(ck~=0) ))*1e-4));
%  Droo=cd;
%  eval(['!md ',Droo,'\',Dire]);
%  Dsav=[Droo,'\',Dire]; disp(' chain: layskip '), keyboard
%end

ichain=1;
if length(find(layskip==1))==length(layskip)
 ichain=0;
end

%' verifica flags', keyboard
%modifica
% if exist('icousav')
%  icoustor=max(icousav);
% else
  icoustor=0;
% end
 icmem=0;
 icaco=1;
 ifplatot=1*isem;
 istr=1;
 T=1;
if ifr==1
% 'entra chain', keyboard
end
if ifp==-11
 'inizio loop strati in chain', pausak
end

if ichain==1
 if ifp==-11
  istr
  pausak
 end
 
 
 for istr=iloo
  
% istr, pausak
  DelT_z0=0;
  icop=fmlsi(istr,1);
  ncop=fmlsi(istr,2);
  %fmlsin=[fmlsin; [icop ncop]];
  
  
  if ifp>0 | ifp==-11
  '   icop  ncop', [icop ncop]
  '   Li  ai  real(ni)  ', [Li(istr) ai(:,istr)' real(ni(:,istr))'], pausak
  end

  istr0=istr
  pausak
  %' set istr0', keyboard
  if itutmir==1
   if ncop~=0
    Ncop=ncop;
   else
    Ncop=1;
   end
  else
   Ncop=1;
  end

  iziv=1;
  ncopF=1;
  istrctv=0;
  ncopco=0;

  ifplaco=1*isem;
  if ifp==-11
   istr
   Ncop
   pausak
  end
  

%  'ncop, icop', [Ncop, icop], keyboard
  
%' Ncot' , keyboard
if igau==4
 zreach=zsopra+sum(zrif([1:istr-1+icop]));
 if zreach<=zedis(1) | zreach>=zedis(end)
  Ncop=1;
  isomv=zreach-sum(Li(istr-1+[1:icop]));
 end
end

if Ncop>1
 %Ncop
 %'Ncop'
 if isto==1, keyboard, end
end
  for istrct=1:Ncop


  
  DelT_z0=0;
   istr=istr0;
if Ncop>1
 istr
 'Ncop istr'
 if isto==1, keyboard, end
 
end   
  %' set istr0 coppie', keyboard
   
   if ncopco==0
    Tc=IdeOon;
    if ifp==-11
     ' Tc = I ', pausak
    end
   end
   ncopco=ncopco+1;
    if ifp==-11
     ncopco
     ' ncopco ', pausak
    end

% qui inizio per struttura ripetuta
%
   for istrc=1:icop
    if Ncop==1 & istrct==1
     fmlsin=[fmlsin; [icop ncop]];
    else 
     fmlsin=[fmlsin; [1 0]];
    end   
    if ifp==-11
     disp(' for: istr, istrc'), [istr istrc], pausak
    end
    DelT_z0=0;
    if igau==4
     isomv=isomv+Li(istr);
     %istr, pausak
     icond=find(isomv>zedis(1) & isomv<zedis(end));
     KTe=0;
     KTez=0;
     izi=0;
     if icond==1
     [du,izi]=min(abs(isomv-zedis(2:length(zedis))));
     DelT_z0=Tdis(end,izi);
     %isomv
     %'Delz ', pausak
%     'Delz ', keyboard
      KTe=reshape(KTemp(:,:,izi),si2);
      if iztm==1
       KTez=reshape(KTempz(:,:,izi),si2);
      else
       KTez=0;
      end
     end
%      if length(find(imag(KTe)~=0))>0
%       DelT_z0=real(DelT_z0)+j*imag(DelT_z0)*pim;
%       KTe=real(KTe)+j*imag(KTe)*pim;
%       KTez=real(KTez)+j*imag(KTez)*pim;
%      end

    end
    
    if ifp==-11
    [istr icaco istrc istrct], pausak
    end
    if icaco==1
     ifiezsav=ifiez;
     ifiez=0;
     if ifp==-11
      Tsave=Tc;
     end
%     [istr icoustor]
%     ' prima di emme', pausak
     eval(emme)

% memorizzo metallo     
    if length(istfie)>0
     if length(find(istr==istfie))==1
      icmem=icmem+1;
      Tmeduf(:,:,icmem)=Oo;
%'memorizzo metallo', pausak
     end
    end     
     
     ifiez=ifiezsav;

% per campo longitudinale

     if ifiez==1
      if icop>1 & Ncop==1
       Tmi{istrc}=Oo;
      end
%      coe_zi
%      'entro coe_zin', keyboard
      if length(Lizi)>1
       coe_zin
      else
       coe_zi
      end
     end  %ifiez

     if istrc==1
      Tcm=Tc;
     end
    end
%     istr
%    if ifast==0
%     mapab(Tc), pausak
%    else
%     'ifast'
%     pausak
%     mapab(reshape(Tstor(:,:,istr),40,40)), pausak
%    end
    if ifp==-11
     ' vedi Tsave ', pausak
    end
%    istr=istr+1;
if Ncop>1
 'Ncop dopo'
 if isto==1, keyboard, end
 
end
   end  %istrc

   if igau==4 & itutmir==1
    if Ncop>1
     if izi~=iziv | istrct==Ncop
      icaco=1;
     else
      icaco=0
      %keyboard
     end
     if icaco==1
      ncopF=ncopco;
      ncopco=0;
     end
    else
     icaco=1;
     ncopco=0;
    end
    iziv=izi;
    istrctv=istrct;
   else
    icaco=1;
    ncopF=ncop;
   end
   if ifp==-11
     ' IMPOR '
     ncop, pausak
   end

   if icaco==1

    if ncop<=1
     if imem==0
      T=Tc*T;
     else
      T=prodmat(Tc,T,ifplatot,Pust);
     end
    else
    if ifp==-11
     ' prima di potensa' , 
       if isto==1, keyboard, end
     
     
    end

     if ncopF-fix(ncopF)==0
      if imem==0
       Tco=Tc^ncopF;
      else
       Tco=powermat(Tc,ncopF,ifplaco,Pust);
      end
      if ifiez==1
%        ' specctio ', keyboard
        for inpa=2:ncopF
         for inpai=1:icop
          pusta=istr-icop-1+inpai;
          dos=Li(pusta);
          shtr=shai(:,pusta);
          nitr=ni(:,pusta);
          Oo=Tmi{inpai};
%          coe_zi
           if length(Lizi)>1
            coe_zin
           else
            coe_zi
           end
         end
%        dis_fz
%        ' dopo paia specctio ', keyboard
        end
      end  %ifiez
      if ifp==-12
       ' dopo specctio ', keyboard
      end
%       ncopF
%      ' ncopF'
%      keyboard
     else
      if imem==0
       Tco=Tcm*Tc^fix(ncopF);
      else
       Pow=powermat(Tc,fix(ncopF),ifplaco,Pust);
       Tco=prodmat(Tcm,Pow,ifplaco,Pust);
      end
%       ' Pow', keyboard
     end

     if imem==0
      T=Tco*T;
     else
      T=prodmat(Tco,T,ifplatot,Pust);
     end
    end
   end  %icaco
%   disp(' istrct '), istrct
%   pausak
  end  %istrct
%Li(istr), pausak

 end   %istr
% clear Tcm Tco Tc
%  disp(' chain crit vecchio'), keyboard

 if ifp==-11 | ifp==-12
  disp(' chain crit'), keyboard
 end

 %if ilaski==0
 if ilaski==10
   if pol==pvet(1)
%    if icoustor~=istr-1

    if exist('icousav')
%     if length(icousav)~=istr-1
     if length(icousav)~=istr-1
      if icoustor~=istr-1
       icoustor=icoustor+1;
       icousav(istr-1)=icoustor;
      end
     else
      icoustor=icousav(istr-1);
     end
     icoustor=icousav(istr-1);
    else
     icoustor=icoustor+1;
     icousav(istr-1)=icoustor;
    end
   else
      icoustor=icousav(istr-1);
   end
%'   qui iprr', keyboard
%if ilaymem(istr-1)==0
%' if ilaymem(istr)==0  cambiato in chain_i_vecchio'
%iprr=1;
%if ilaymem(istr-1)==0
% iprr=0;
%end

% vecchio
if ilaymem(istr-1)==0
 iprr=0;
else
 iprr=1;
end


%iprr=0;

 if iprr==0
  if iTsav==0
   if length(T)>1
     if ifp~=-4 & ick==1
      istr
      'salvo T', keyboard
     end
    Tstor(:,:,icoustor)=T;
   end
  else
   Tstof=T;
   if ispeed==1
    if length(T)>1
     eval(['save ',Dsav,'\', nTstof, num2str(icoustor),' Tstof']);
    end
   end
  end
 end  %iprr

 end

 if ifp==-11
  disp(' chain dopo'), keyboard
 end
end  %ichain

%ifp=ifpss;

  %disp(' chain dopo'), keyboard

 Tdu=IdeOon;
  if ick==1
' % Tdu=1; ', keyboard
  end

if length(Li)==0
 Tb=1;
 jsau=0;

else 

 if iTsav==0
  sj=size(Tstor);
  if length(sj)==3
   jsau=sj(3);
  else
   jsau=1;
  end
 else
  jsau=icoustor;
 end

end

 if ifp==-12
  disp(' jsau '), keyboard
 end

if ifp~=-4 & ick==1
' icousav', keyboard
end

%' sono in icrit: chain_vecchio', keyboard
 ficri=find(icrit==1);
 if isfield(Ps,'igacrit')
  if Ps.igacrit==0
  'RESETTO ICRIT'
   ficri=[];
   icrit=icrit*0;
  end
 end
 if dovesono==1
  ficrit=ficri;
 else
  ficrib=ficri;
 end 
 

%ifiez, keyboard

if length(find(icrit>0))>0 & ifiez~=1

%if exist('icrit')==1
 
 %sha_all=shavet(2:end-1); 
 %clear nimi
 %for ks=1:length(nitot)
 % nim=find(nitot(ks,:)>0);
 % nimi(ks)=min(nitot(ks,nim));
 %end 
 %ficr=find(sha_all>3 & nim<2 & aitot>0);
 %jsa=ficr(end);
 %icrit(jsa+1:end)=0;
% ' qui cont nuovo prima', keyboard
if dovesono==1
  jsa=ficri(end); 
  jsau=icoustor; 
else
  jsa=0; 
  jsau=ficri(1)-1; 
end
 % [Gacritn,Trcrit]=Gam_critU(Tstor,Ga1,Mcrit,ficri);
 if ifp==-10
%  ' qui cont nuovo', keyboard
 end
 if dovesono==1
  Ga_old=Ga1;
  Ga_inf=Ga1;
 else
  Ga_old=Ga2;
  Ga_inf=Ga2; 
 end  

%'Gacrit', keyboard

% [Gacrit,Trcrit,Trc,Grc]=Gam_critUn(Tstor,Ga1,Mcrit,ficri,fmlsi);
ficris=ficri;
if igraef_new==2
 fisha=find(shai==6);
 ficri(fisha)=-ficri(fisha);
end 


[Gacrit,Trcrit,Trc,Grc]=Gam_critScattHCG(Ga_inf,Mcrit,ficri,fmlsi,dovesono,icpo);


% [Gacrit,Trcrit]=Gam_critU(Tstor,Ga1,Mcrit,ficri);
 if ifp==-10
%  ' sono in icrit: chain_vecchio DOPO', keyboard
 end
 if dovesono==1
  Trcritu{icpo}=Trcrit;
  Ga1=Gacrit;
 else
  Ga2=Gacrit;
  Trcritb=Trcrit;
 end  
% ' sono in icrit: chain_vecchio', keyboard
else 
 jsa=0;
end
' jsau Mcri ', keyboard
return

if jsau>0
% jsa=0;
 while jsa<jsau

  jsa=jsa+1
  nrig=find(icousav==jsa);
  if length(nrig)>1
%   ' nrig in Chain_i ', keyboard
  end
  nstrat=fmlsi(nrig,1);
  if ifp==-11
  disp(' in while 1')
  jsa
  nstrat
  nrig
  pausak
  end
  if ilaymem(nrig)==0
   nstrat=1;
  end
  if nstrat==1
   if iTsav==0
    Tdu=Tstor(:,:,jsa)*Tdu;
   else
%    eval([' load ',nTstof,num2str(jsa)]);
    if ispeed==1
     eval([' load ',Dsav,'\',nTstof,num2str(jsa)]);
    end
    Tdu=Tstof*Tdu;
   end
  else
   nmirro=max(fmlsi(nrig,2))
%   keyboard
   if ifp==-11
    'nmirro'
    nmirro
    pausak
    keyboard
   end
   Tmirro=IdeOon;
   for kmir=1:nstrat
    if iTsav==0
     Tdum=Tstor(:,:,jsa);
    else
     if ispeed==1
      eval([' load ',Dsav,'\',nTstof,num2str(jsa)]);
     end
     Tdum=Tstof;
    end
    Tmirro=Tdum*Tmirro;
    jsa=jsa+1;
    if ifp==-11
     disp(' in while 2')
     jsa
    end
   end

%' nmirro', keyboard
   Pow=Tmirro^nmirro;
   Tdu=Pow*Tdu;
    if ifp==-11
    'Tmirro ', keyboard
     disp(' in while 3')
     jsa
     pausak
    end
   jsa=jsa-1;
  end
%  istr
%  pausak
%  if length(find(istfie==istr))==1
%   icmem=icmem+1;
%   Tmeduf(:,:,icmem)=T;
%   disp(' memorizzo  Tmef')
%   keyboard
%  end
%  jsa
%  pausak
%  if length(find(istfie==jsa))==1
  
%  puf=find(ilaymem==1);
%  if length(puf)==0
%   puf=-1;
%  end
%  if length(istfie)>0
%  if puf==istfie
%   icmem=icmem+1;
%   Tmeduf(:,:,icmem)=Tdu;
%   if ifp~=-4 & ick==1 , disp(' memorizzo  Tmef'), keyboard, end
%  end
 end

else
 if length(Li)>0
 if iTsav==0
  Tdu=Tstor;
 else
  Tdu=Tstof;
 end
 end
end

if ifp==-12
  disp('fine  Tstor '), keyboard
end
%  disp('fine  chain Tstor '), keyboard
% clear T Mo
