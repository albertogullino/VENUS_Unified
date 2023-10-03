ideb=0;
%' chain', keyboard
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
 fisem=find(ai(1,:)~=0);
% asem=zeros(1,length(ai));
 asem=zeros(1,length(ai(1,:)));
 asem(fisem)=ai(1,fisem);
% ilaymem=find(diff([asem 0])<0 & shailoc>1);
% ilaymem=find(diff([asem 0])<0 );
% ilave=find(asem~=0 & shailoc>1);
% ilave=find(asem~=0 & shailoc>0 & Li'~=0);
 ilave=find(asem~=0 & shailoc>0 );
 ilaymem=zeros(1,length(ai));
 ilaymem(ilave)=1;
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
  ilave=find(asem~=0 & shailoc>0);
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
% 'entra chain', keyboard

if ichain==1
 while istr<=length(Li)
  if ideb==1
  if istr==23
   ' fine specchio ', keyboard
  end
  end
  DelT_z0=0;
  Dos=Li(istr);
  icop=fmlsi(istr,1);
  ncop=fmlsi(istr,2);
  if ifp>0
  '   icop  ncop', [icop ncop]
  '   Li  ai  real(ni)  ', [Li(istr) ai(:,istr)' real(ni(:,istr))'], pausak
  end

  istr0=istr;
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

  for istrct=1:Ncop

  DelT_z0=0;
   istr=istr0;
   if ncopco==0
    Tc=IdeOon;
   end
   ncopco=ncopco+1;

   for istrc=1:icop
    if ifp>=1
     disp(' for: istr, istrc'), [istr istrc], pausak
    end
%     disp(' for: istr, istrc'), [istr istrc], pausak
    DelT_z0=0;
    if igau==4 & sT>0
     isomv=isomv+Li(istr);
%     pausak
%     [du,izi]=min(abs(isomv-zedis(2:length(zedis))));
     [du,izi]=min(abs(isomv-zedis));
     DelT_z0=Tdis(end,izi);
%     'Delz ', keyboard
      KTe=reshape(KTemp(:,:,izi),si2);
      if iztm==1
       KTez=reshape(KTempz(:,:,izi),si2);
      else
       KTez=0;
      end

    end
    if ifp==-11
    [istr icaco istrc istrct], pausak
    end
    if icaco==1
     ifiezsav=ifiez;
     ifiez=0;
     eval(emme)
     ifiez=ifiezsav;

% per campo longitudinale

        if ifiez==1
         if icop>1 & Ncop==1
          Tmi{istrc}=Oo;
         end
%         coe_zi
%         'entro coe_zin', keyboard
         if length(Lizi)>1
          coe_zin
         else
          coe_zi
         end
        end  %ifiez
% fine per campo longitudinale

     if istrc==1
      Tcm=Tc;
     end
    end

%  'istrc fine ', istr
%  pausak

    istr=istr+1;

   end  %istrc
   if ideb==1
   'fine istrc', pausak
   end

   if igau==4 & itutmir==1
    if Ncop>1
     if izi~=iziv | istrct==Ncop
      icaco=1;
     else
      icaco=0;
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
   if ifp>=1
     ' Ncop = ', ncopF, pausak
   end

   if icaco==1

    if ncop<=1
      if ideb==1
       ' prodotto T', pausak
      end
     if imem==0
      T=Tc*T;
     else
      T=prodmat(Tc,T,ifplatot,Pust);
     end
    else

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
%       dis_fz
%        ' dopo paia specctio ', keyboard
        end
     end  %ifiez
%      ' dopo specctio ', keyboard
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
     end

     if imem==0
      T=Tco*T;
     else
      T=prodmat(Tco,T,ifplatot,Pust);
     end
      if ideb==1
       ' Pow', keyboard
      end
    end
   end  %icaco
   if ideb==1
    disp(' istrct '), istrct
    pausak

    if istrct==23
    'fine specchio ', keyboard
    end
   end
  end  %istrct

 end   %istr
% clear Tcm Tco Tc

 if ifp==-11
  disp(' chain crit'), keyboard
 end

 if ilaski==0
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
if ilaymem(istr-1)==0
 iprr=0;
else
 iprr=1;
end
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

      if ideb==1
       '   qui iprr', keyboard
      end
 end

 if ifp==-11
  disp(' chain dopo'), keyboard
 end
end  %ichain

 Tdu=IdeOon;
  if ick==1
' % Tdu=1; ', keyboard
  end

%'1ui ver', keyboard
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


 if ifp==-11
  disp(' jsau '), keyboard
 end

if ifp~=-4 & ick==1
' icousav', keyboard
end

if jsau>0
 jsa=0;
 while jsa<jsau

  jsa=jsa+1;
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

   Pow=Tmirro^nmirro;
   Tdu=Pow*Tdu;
%   'Tmirro ', keyboard
    if ifp==-11
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
  puf=find(ilaymem==1);
  if length(puf)==0
   puf=-1;
  end
  if puf==istfie
   icmem=icmem+1;
   Tmeduf(:,:,icmem)=Tdu;
   if ifp~=-4 & ick==1 , disp(' memorizzo  Tmef'), keyboard, end
  end

 end
else
 if iTsav==0
  Tdu=Tstor;
 else
  Tdu=Tstof;
 end
end

if ifp==-11
  disp('fine  Tstor '), keyboard
end

% clear T Mo
