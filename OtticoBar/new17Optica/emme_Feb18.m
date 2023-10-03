if ~exist('ioxi')
 ioxi=0;
end 
%' emme_ult ', keyboard
%kv=KKt(Pusas);
ick=1;
ick=0;
if ick==1
if ifp~=-4
istr
if istr==42
 ifpsap=ifp;
 ifp=3;
else
 if exist('ifpsap')
 ifp=ifpsap;
 end
end
end
end  %ick
%ifp=-10
%Kan_gr=0;
ibrea=0;

%    if iLP==1 & istr==5
%     'emmemix 2 LP', keyboard
%     ibrea=1;
%    end
%%%%%%%%%%%%%%%/pim% All'interno di un ciclo su ciascuno strato istr,

         if ifiez==0
          ilayfastistr=ilayfast(istr);
          layskipistr=layskip(istr);
          dos=Li(istr);
         elseif ifiez==1
          dos=Liz;
          ilayfastistr=0;
          layskipistr=0;
         elseif ifiez==2
          ilayfastistr=0;
          layskipistr=0;
         end
if istr==1
%         'mix'
%keyboard
end
ifpsal=ifp;

if shai(1,istr)==0 | Li(istr)==0 | ai(istr)>=0
% ifp=-4;
end
if ioxi==1
  '[ifiez dos ilayfastistr layskipistr]'
  [ifiez dos ilayfastistr layskipistr]
  ifp=5;
  pausak


end

%if si2(1)>200
%' % calcolo M , si2'
%keyboard
%end
global PUrid

%tic
%istr
%dpes=diag(pes);
ilin=0;    % metodo raffinato : linearizza dal punto precedente
%ilin=1;    % metodo + rude: linearizza sempre dal primo punto

if ilin==0
 pufreq=ifr-1;
else
 pufreq=1;
end
PUrid=0;

if length(find(shavet==6))>0
 iproga=1;
else
 iproga=0;
end

 iproga=0;

[Tfas,Ffas]=eltime(Tfas,Ffas,istr);

if ifp>=2
disp(' entro emmestru: istr ')
istr
end

%keyboard
%disp(' entro emmestru: istr ')
%pausak

if istr>1 & ifiez==0
  if layskip(istr-1)==0 & ilaymem(istr-1)==0 & ilaymem(istr)==1
   if pol==pvet(1)
    icoustor=icoustor+1;
    icousav(istr-1)=icoustor;
   else
    icoustor=icousav(istr-1);
   end
%   'inizio'
%  [icustor istr]
%   pausak
   if iTsav==0
    if length(T)==1
     Tstor(:,:,icoustor)=diag(ones(size(Pust))) ;
    else
     Tstor(:,:,icoustor)=T;
    end
   else
    if length(T)==1
     Tstof=diag(ones(size(Pust))) ;
    else
     Tstof=T;
    end
    if ispeed==1
     eval(['save ',Dsav,'\', nTstof, num2str(icoustor),' Tstof']);
    end
%    'salvato'
%    keyboard

   end
   ifplatot=1*isem;
   T=IdeOon;
  end
end



%if layskip(istr)==0
if layskipistr==0


%  iama=iaccvet{istr};
  aitr=ai(:,istr);
  bitr=bi(:,istr);
  pitr=pai(:,istr);
  shtr=shai(:,istr);

  if ity==1
   tytr=tyari(:,istr);
  end

  nitr=ni(:,istr);
  aniloc=ani(istr);
  aniloc1=ani1(istr);
  if ~exist('anir')
   anir=zeros(size(ani));
  end
  ani_gr=anir(istr);  
  if (ani_gr)~=0
   rintz=rint+ani_gr/2;  % eps_z = eps_x
    if ick==1
     'ani in emme', keyboard
    end 
  end
  if ~exist('Kan_gr') 
   Kan_gr=0;
  end


    if aitr==0 & igau~=4
     if imem==1
      ifpla=1*isem;
      if (aniloc~=0 & iany==1 & ianys<=1) | (aniloc1~=0 & ianys==1 & iany<=1)
       ifpla=2*isem;
      end
      if ani_gr~=0
       ifpla=2*isem;
      end
     else
      ifpla=0;
     end
    else
     ifpla=0;
    end
    ifplatot=ifplatot*ifpla;
     ifplaco=ifplaco*ifpla;


%    disp('  istr, ifpla, ifplaco, ifplatot ')
%   [istr, ifpla, ifplaco, ifplatot],
%   pausak

  rint=nitr(1);


% parte per profilo di perdite trasversale


      pim=imag(rint^2);

 DelT_z=0;
 if igau==4

%      if pim~=0
%       'pim diverso 0', keyboard
%      end
  DelT_z=DelT_z0;
  KTep=dpes1*KTe*dpes2/rr^2;
  if iztm==1
   KTepz=dpes1*KTez*dpes2/rr^2;
  else
   KTepz=0;
  end

 else
  KTep=0;
  KTepz=0;
 end



  fai=find(aitr~=0);
  if length(fai)==0
   aidu=0;
   rext=0;
  else
   aidu=aitr(fai(1));
   sgai=1;
   if exist('iGRAD')
    if PAlgrad.Isha==0 | PAlgrad.Iapp==1
     if iGRAD(istr+1)==1
      aidu=abs(aidu);
      sgai=-1;
     end
    end
   end
%   ' aidu', keyboard
   ini=find(aitr==0);
   if length(ini)>0
    riv=nitr(1:ini(1));
   else
    riv=nitr;
   end
%   ' riv vedo ',
%   istr
%   fai'
%   pausak
    rext=riv(end);
    if rext~=1
     rext=rext+DelT_z/(2*rext);
    end
  end
%  if istr==32
%  ' istr=32', keyboard
%  end
  if ibrea==1
   ' qui' , keyboard
  end
%  if aidu==0
  if aidu<0
   dos=0;
  end
%%%%%db%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Se lo strato presenta una struttura trasverale non banale
% if aidu~=0 & Li(istr)~=0

% if aidu~=0 & dos~=0
 if aidu>0 & dos~=0
     if ifp>=2
        display('sono entrata della parte con struttura trasversa'), pausak
     end

%    dos=Li(istr);
    Delta=ani(istr)/rr^2;
    Delta1=ani1(istr)/rr^2;
    Delta_gr=anir(istr)/rr^2;
    if abs(Delta_gr)>0 
     if ick==1
      'emmemix', keyboard
     end
    end
    if ifp>0
       disp('[spessore nint next]')
       [dos rint rext]
       disp('[a n]')
       [[aitr(fai); 0] riv]
       pausak
    end

%%%%%%%%% a)

        Kost=0;
        Kostz=0;
        for ianu=1:length(fai)   %(ciclo su ogni sezione trasversale omogenea)
            at0=aitr(ianu);
            bt0=bitr(ianu);
            pt0=pitr(ianu);
            sh0=abs(shtr(ianu));
            if sh0==7
             sh01=1;
            else
             sh01=sh0;
            end


%if sh0==1
%ifp=3
%else
%ifp=-10
%end
%            pausak


            if sh01>=5
             if ity==1
              ty0=tytr(ianu);
             end
             depepf=(riv(ianu)^2-riv(end)^2)/rr^2;
             depepfz=-rr^2/riv(ianu)^2+rr^2/riv(end)^2;
%             depepf=riv(ianu)^2/rr^2-1;
%             depepfz=1-rr^2/riv(ianu)^2;
            else
             depepf=(riv(ianu)^2-riv(ianu+1)^2)/rr^2;
             depepfz=-rr^2/riv(ianu)^2+rr^2/riv(ianu+1)^2;
            end

            ipri=0;
            if sh01==1
             ipri=find(aytot(:,sh0)==at0);
%             if length(ipri)>1
%              ipri=ipri(1);
%             end
            elseif sh0>1 & sh0<4
             ipri=find(aytot(:,sh0)==at0 & axtot(:,sh0)==bt0);
            elseif (sh0==4 | sh0==6)
             ipri=find(aytot(:,sh0)==at0 & axtot(:,sh0)==bt0 & pdtot(:,sh0)==pt0);
%             if length(ipri)>1
%              ipri=ipri(2);
%             end
            elseif sh0==5
             ipri=find(aytot(:,sh0)==at0 & tytot==ty0);
            elseif sh0==8
             ipri=1;
            end
%           ' selezione e_mix', keyboard
            Kr=0;
%            istr
%            ' wmme mix', pausak
            if ipri~=0
             iama=Iacc{sh0};
             PUrid=PUrD{sh0};
%             iama=[];
%             pausak
             KOS=Kos{sh0};


            if igau==5 
             KTep=dpes1*KTemp*dpes2/rr^2;
	      if iztm==1
	       KTepz=dpes1*KTempz*dpes2/rr^2;
	      else
	       KTepz=0;
             end
            end

             if  length(size(KOS))==3
%              Kr=dpes1*reshape(KOS(:,:,ipri),si2)*dpes2*depepf;
              Kr=dpes1*KOS(:,:,ipri)*dpes2*depepf;
             else
%              Kr=KOS*diag(pes)*depepf;
              Kr=dpes1*KOS*dpes2*depepf;
             end
%            disp(' emme_mix ox')
%            keyboard
%             clear KOS
            else
             Kr=0;
            end

            if ifp>=1
             ' ferma: at0 bt0 pt0 sh0 ipri = ',[at0 bt0 pt0 sh0 ipri], keyboard
            end

            Kost=Kost+Kr;
            if iztm==1
              KOS=Kosz{sh0};
             iama=Iacc{sh0};
             PUrid=PUrD{sh0};
%             iama=[];
%             pausak
              if  length(size(KOS))==3
%               Krz=dpes1*reshape(KOS(:,:,ipri),si2)*dpes2*depepfz;
               Krz=dpes1*KOS(:,:,ipri)*dpes2*depepfz;
              else
%               Krz=KOS*diag(pes)*depepfz;
               Krz=dpes1*KOS*dpes2*depepfz;
              end
              clear KOS

              Kostz=Kostz+Krz;
            end
        end  %(fine ciclo sulle sezioni trasversali)
            if ifp>=2
%            if istr==2
              disp('[ianu depepf]')
              [ianu depepf]
              figure, surf(abs(Kost)), shading('interp'), view(0,90),
              title(' Kost ')
              colorbar, pausak
              disp(' emmestru ox circ')
              pausak
            end

        if ifp>=20
           figure,
           surf(abs(Kost)), shading('interp'), view(0,90),
           title(' KOST ')
           colorbar, pausak
        end



         depepo=(rext^2-rr^2)/rr^2;
         depepoz=1-rr^2/rext^2;
        if at0>ra0
         depepop=(rint^2-rr^2)/rr^2;
         depepozp=1-rr^2/rintz^2;
         IDE=Madd0*depepop+(Idelta-Madd0)*depepo;
         if iztm==1
          IDEz=Maddz0*depepozp+(Ideltaz-Maddz0)*depepoz;
         end
%         ' emme mesa ', keyboard
        else
         Madd=0;
         Maddz=0;
         IDE=Idelta*depepo;
         IDEz=Ideltaz*depepoz;
        end
        ck=-j*kcav*dos;
        Belo=-j*be*kcav*dos;
        if iztm==1
%           KOt=ck*(Idelta*depepo+Delta*Kan+Delta_gr*Kan_gr+Delta1*Kan1+Kost+KTep);
           KOt=ck*(IDE+Delta*Kan+Delta_gr*Kan_gr+Delta1*Kan1+Kost+KTep);
%           if Delta_gr~=0
%           KOz=0;
%           ' Delta_gr'
%           KOz
%           else
%           KOz=ck*(Ideltaz*depepoz+Kostz+KTepz);
           KOz=ck*(IDEz+Kostz+KTepz);
%           KOt(1:20,:)=0; 
%           KOt(:,1:20)=0; 
%           KOz(1:20,:)=0; 
%           KOz(:,1:20)=0; 
%           end

%           Ksu=KOt+KOz;
%           Kdi=KOt-KOz;
%           P=[Ksu Kdi; -Kdi -Ksu];    %(matrice totale con ff fb bf bb)
%            P=[KOt+KOz  KOt-KOz; -(KOt+KOz)  -(KOt-KOz)];
            P=[KOt+KOz  KOt-KOz; -(KOt-KOz)  -(KOt+KOz)];
% ' PPP ', keyboard
        else
%           KOt=ck*(Idelta*depepo+Delta_gr*Kan_gr+Delta*Kan+Delta1*Kan1+Kost+KTep);
%           KOt1=ck*(IDE+Kost);
           KOt=ck*(IDE+Delta_gr*Kan_gr+Delta*Kan+Delta1*Kan1+Kost+KTep);
           P=[KOt KOt; -KOt -KOt];
        end
       if ifp==1
        [istr; aitr], 'emme_mix', pausak
       end

%        Mo=diag([Belo; -Belo])+P;
        P=diag([Belo; -Belo])+P;
%        if nitr(2)==1.6
%         'ver exp P '
%         keyboard
%        end
%    if abs(Delta_gr)>0
%     'emmemix 2 vet dentro', keyboard
%    end
%    if iLP==1 & istr==5
%     'emmemix 2 LP dentro', keyboard
%    end

%          if ilayfast(istr)==0
          if ilayfastistr==0

%            ifpla
%            if ifpla==2
%             ifpla=0;
%            end
%           Oo=expm_mio(Mo,ifpla,iama,Pust);
%$           if length(inuo_bas)>0

%    if abs(Delta_gr)>0
%     'Oo', keyboard
%    end
%           if sh0==1
%            ' verifica ', keyboard
%           end
           if inuo_bas==0
%            size(Oo)
            expm_miu(P,ifpla,iama,Pust);
%            size(Oo)
%            ' emme ver', pausak
           else
%            ' SET Oo 1', keyboard
            Oo=expm_mio(P,ifpla,iama,Pust,PUrid);
           end
           if ifp>0 & ick==1
           'Oo', keyboard
           'Oo', keyboard
           'Oo', keyboard
           'Oo', keyboard
           end

           freqs=freq;
%           if ifr==2
%            'qui non sem', pausak
%        Oos=Oo; save sa Oos
%        load sa, mapab(Oos-Oo)
%           end

          else  %ilayfast
%           'qui', pausak
           icousa=icousav(istr);
%           if ifr==2
%            'qui non sem', pausak
%            keyboard
%           end
           if iTsav==0
%            Oo=Tstof{icousa}+Tstof{icousa}*Mo(Pust,Pust)*(freq-fre_camp(pufreq));
%            ' SET Oo 2', keyboard
            Oo=Tstof{icousa}+Tstof{icousa}*P(Pust,Pust)*(freq-fre_camp(pufreq));
           else
            if ispeed==1
             eval([' load ',Dsav,'\',nTstof,num2str(icousa)]);
            end
%            Oo=Tstof+Tstof*Mo(Pust,Pust)*(freq-fre_camp(pufreq));
%            ' SET Oo 3', keyboard
            Oo=Tstof+Tstof*P(Pust,Pust)*(freq-fre_camp(pufreq));
           end

          end  %ilayfast





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Se la struttura trasversa e' banale   (omogenea)

 else % aitr  ai

   if ifp>=2
      disp('sono entrata della parte senza struttura trasversa'), pausak
   end
%      disp('sono entrata della parte senza struttura trasversa'), pausak

%   dos=Li(istr);
   Delta=ani(istr)/rr^2;
   Delta1=ani1(istr)/rr^2;
   Delta_gr=anir(istr)/rr^2;
     if ifp>=2
       disp('[spessore nint larghezza anisotropia]')
       [dos rint ai(istr) Delta]
       pausak
     end
     if rint~=1
      rintv=rint;
      rint=rint+DelT_z/(2*rint);
     end
     rintz=rint;

    if dos~=0
       depe=(rint^2-rr^2)/rr^2;
%       depe, keyboard
       depez=1-rr^2/rintz^2;
       ck=-j*kcav*dos;
       Beli=ck*be;
       Mod=([Beli; -Beli]);
       if iztm==1
          KOt=ck*(Idelta*depe+Delta*Kan+Delta_gr*Kan_gr+Delta1*Kan1+KTep);
%           if Delta_gr~=0
%           KOz=0;
%           ' Delta_gr'
%           KOz
%           else
           KOz=ck*(Ideltaz*depez+KTepz);
%           end

%          Ksu=KOt+KOz;
%          Kdi=KOt-KOz;
%          P=[Ksu Kdi; -Kdi -Ksu];
%          P=[KOt+KOz  KOt-KOz; -(KOt+KOz)  -(KOt-KOz)];
          P=[KOt+KOz  KOt-KOz; -(KOt-KOz)  -(KOt+KOz)];
       else
          KOt=ck*(Idelta*depe+Delta*Kan+Delta_gr*Kan_gr+Delta1*Kan1+KTep);
          P=[KOt KOt; -KOt -KOt];
       end
          Pk_sav=P;
    if abs(Delta_gr)>0
     if ifr==1 & ifp~=-4 & ick==1
      pr=Idelta*depe+Delta_gr*Kan_gr;
     'emmemix 1', keyboard
     end
    end

       if ifp>=2
         display('matrice K'), pausak
         figure,  surf(abs(KOt)), shading('interp'), view(0,90),
              title(' KOt ')
         colorbar, pausak
       end

%            Mo=diag(Mod)+P;
            P=diag(Mod)+P;
%            P=diag(Mod)+P(Pust);
%'expm_mio 2', keyboard
%pausak

      %      ifpla
      %      if ifpla==2
      %       ifpla=0;
      %      end
          if ilayfastistr==0
%          ilayfastistr=ilayfast(istr);
%          if ilayfast(istr)==0

           iama=[];
%           Oo=expm_mio(P,ifpla,iama,Pust);
%           if length(inuo_bas)>0
           if inuo_bas==0
%            size(Oo)
%            ' SET Oo 4',
%            if ifr==2
%            keyboard
%            end

            expm_miu(P,ifpla,iama,Pust);
%            size(Oo)
%            ' emme ver', pausak
           else
%            ' SET Oo 5', keyboard
            Oo=expm_mio(P,ifpla,iama,Pust,PUrid);
%            ' SET Oo 5 dopo', keyboard
           end

%'vedo Oo', keyboard


           freqs=freq;

          else  %ilayfast

           icousa=icousav(istr);
           if iTsav==0
%            Oo=Tstof{icousa}+Tstof{icousa}*Mo(Pust,Pust)*(freq-fre_camp(pufreq));
%            ' SET Oo 6', keyboard
            Oo=Tstof{icousa}+Tstof{icousa}*P(Pust,Pust)*(freq-fre_camp(pufreq));
           else
            if ispeed==1
             eval([' load ',Dsav,'\',nTstof,num2str(ficousa)]);
            end
%            Oo=Tstof+Tstof*Mo(Pust,Pust)*(freq-fre_camp(pufreq));
%            ' SET Oo 7', keyboard
            Oo=Tstof+Tstof*P(Pust,Pust)*(freq-fre_camp(pufreq));
           end

          end  %ilayfast
    if ifp==-11
     ' fine strato piano '
     keyboard
    end

   else
%            ' SET Oo 8', keyboard
    Oo=1;
   end

  end  %ai
%     istr
%     figure, plot(diag(abs(Oo))), pausak
if ifp>0 | ifp==-11
'istr=', istr
' qui emme_ult',% keyboard
pausak
end

     if imem==0
      Tc=Oo*Tc;
     else
      Tc=prodmat(Oo,Tc,ifpla,Pust);
     end
if ifp>0 | ifp==-11
'istr=', istr
' qui dopo emme_ult',% keyboard
pausak
end

  if ifp>=2 & length(Oo)>5
   figure, surf(abs(Oo)), shading('interp'), view(0,90),
   title(' Oo ')
   colorbar, pausak
  end

else
%            ' SET Oo 9', keyboard
  Oo=1;

end  %layskip

if ifp~=-4 & ick==1
[istr]
'emme', pausak
end

if ifiez==0
if layskip(istr)==0 & ilaymem(istr)==1 & (ifr==1 | ilin==0)
   if pol==pvet(1)
    icoustor=icoustor+1;
    icousav(istr)=icoustor;
   else
    icoustor=icousav(istr);
   end
%   'fine'
%   [icoustor istr]
%   pausak

   if iTsav==0
    if length(Tc)==1
     Tstor(:,:,icoustor)=diag(ones(size(Pust)));
    else
     if ncop==0
      Tstor(:,:,icoustor)=Tc;
     else
%      ' Tstor ', pausak
      Tstor(:,:,icoustor)=Oo;
     end
    end
   else
    if length(Tc)==1
     Tstof=diag(ones(size(Pust)));
    else
     if ncop==0
      Tstof=Tc;
     else
      Tstof=Oo;
     end
    end
    if ispeed==1
%     ['save ',Dsav,'\', nTstof,num2str(icoustor),' Tstof']
     eval(['save ',Dsav,'\', nTstof,num2str(icoustor),' Tstof']);
    end
%    'salvato 1'
%    keyboard

   end
   prga_ret

%   [icoustor istr]
%   disp(' emme_mix fin: Tstor'), keyboard
   ifplatot=1*isem;
%   T=1;
%   Tc=1;
   T=IdeOon;
   if ifp<=-11
     [istr, istrc]
     ' acktu!!!!!! metto Tc=1 ', pausak
   end
   Tc=IdeOon;
end
end  %ifiez

[Tfas,Ffas]=eltime(Tfas,Ffas,istr);

ifp=ifpsal;
%istr
%' % dopo calcolo M ', keyboard
if ifp>0
' % dopo calcolo M '
 if ifp>2
  keyboard
 end
end
