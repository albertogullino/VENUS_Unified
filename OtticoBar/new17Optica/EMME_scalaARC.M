iautm=1;
           ianisotropia=1;
icriti=icrit(istr);
istocont=0;
%'inizion emme', keyboard
if istr==1700
' entro emme scale', keyboard
end
%kcav=1.074717021057574e+001
global Dsav
Oo=IdeOon;
iSET=0;
%if istr==88
% iSET=1;
%end
if iSET==1
'inizion emme', keyboard
end
sh0=-1000;

if pol==-1
 Idelta=Idelta_m;
 Ideltaz=Ideltaz_m;
else
 Idelta=Idelta_p;
 Ideltaz=Ideltaz_p;
end
%iSET0=0;

if length(ilaymem)==0
  ilaymem(istr)=0;
end

if ~exist('igraef_new')
  if isfield(Ps,'igraef_new')==1
    igraef_new=Ps.igraef_new;
 else
   igraef_new=0;
 end 
end

iSET0=0;
if igraef_new==-1
 iSET0=1;   % si ferma
end 
%iSET0=0;
%save emme
%'salvato emme', keyboard
if exist('Par')
if length(Par)==2
 dpes1=Par{2}.dpes1;
 dpes2=Par{2}.dpes2;
 Idelta=Par{2}.Idelta;
 Ideltaz=Par{2}.Ideltaz;
 Madd0=Par{2}.Madd0;
 Madd0z=Par{2}.Maddz0;
 be=Par{2}.be;
end 
end

iredot=0;
icontrolla=0;
dialef=1;
istr_ret_curv=0;
igra=0;
if ~exist('ioxi')
 ioxi=0;
end 
if dovesono==2
% istr
% ' emme_any0 ', keyboard
end
%kv=KKt(Pusas);
ick=1;
ick=0;
ickset=0;
if ick==1
if ifp~=-4
istr
if istr==420
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
          sh0=shai(1,istr);
          if ~exist('SPP')
           SPP=0;
           nScalini=1;
          end
          if sh0==11 & SPP==1 & Li(istr)~=0
%           'layfast', keyboard
		if ScalaARC==1
		 ScaInc=nScalini/2;
		else
		 ScaInc=nScalini;
		end 
           ianisotropia=0
           dos=Li(istr)/ScaInc; 
           if shai(1,istr+1)==11
            NsMet=Par.NsMet;
            NsNomet=nScalini-NsMet;
            dos=Li(istr)/NsNomet;
           else
            if isfield(Par,'NsMet')
             if Par.NsMet>0
              dos=Li(istr)/NsMet;             
             end
            end 
	   end 
           ilayfastistr=0;
           layskipistr=0;
          end
          

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

%ilayfastistr, pausak


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
   'inizio'
  [icustor istr]
   pausak
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
  prtr=prosi(:,istr);
%  'prtr', keyboard
  aitr=ai(:,istr);
  fier=find(aitr>0);
  aierr=aitr(fier);
  if length(find(diff(aierr)<0))>0
   ' in Emme_any0_11'
   ' File errato !!!!  raggi in ordine errato ', keyboard
  end
  bitr=bi(:,istr);
  pitr=pai(:,istr);
  shtr=shai(:,istr);
  


  if ity==1
   tytr=tyari(:,istr);
  end

  nitr=ni(:,istr);
  
  if length(find(imag(nitr)<-1))
%   ifp=3
%   keyboard
  end
  aniloc=ani(istr);
  aniloc1=ani1(istr);
  if ~exist('anir')
   anir=zeros(size(ani));
  end
  ani_gr=anir(istr);  
% iSET=1; 
 if ani_gr~=0
   iSET=iSET0;
%  'retirnoc', keyboard
 end
  if ~exist('Kan_gr') 
   Kan_gr=0;
  end
  if ~exist('Kanz_gr') 
   Kanz_gr=0;
  end

    if aitr==0 & igau~=4 & prtr>0
%    if aitr==0 & igau~=4 & prtr>0
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

%'qui ipla', keyboard

%    disp('  istr, ifpla, ifplaco, ifplatot ')
%   [istr, ifpla, ifplaco, ifplatot],
%   pausak

  fai=find(aitr~=0)';
  if length(fai)==0
   fai=1;
  end 

  fi6=find(shtr==6);
  if length(fi6)==1 & igraef_new~=1 
   fai=fi6;
  end

  rint=nitr(fai(1));
  rintz=rint;

%' parte per profilo di perdite trasversale', keyboard


      pim=imag(rint^2);

 DelT_z=0;
 if igau==4 
  KmeT=mean(mean(KTe));

 else
  KmeT=0;
 end 
 
 if ifr==1
 istr
 KmeT
 %pausak
 end
 
 if KmeT~=0
  if ~exist('FiMAT')
   subsem
  end
      if length(find(imag(KTe)~=0))>0
       DelT_z=real(DelT_z0)+j*imag(DelT_z0)*pim;
       KTe=real(KTe)+j*imag(KTe)*pim;
       KTez=real(KTez)+j*imag(KTez)*pim;
      else
       DelT_z=real(DelT_z0);
      end
%      if pim~=0
%       'pim diverso 0', keyboard
%      end

  KTep=dpes1*KTe*dpes2/rr^2;
  if iztm==1
   KTepz=dpes1*KTez*dpes2/rr^2;
  else
   KTepz=0;
  end
%[  istr, isomv, dos, izi]
%  'KTE', pausak
 else
  KTep=0;
  KTepz=0;
 end




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
   
   ini=find(aitr~=0)';
   if length(ini)>0
    riv=[nitr(ini); nitr(ini(end)+1)];
   else
    riv=nitr;
   end   
    riv=nitr;   
   %' riv vedo ', keyboard
%   istr
%   fai'
%   pausak
  if length(riv)==1
    rext=riv(end);
  else
    rext=riv(fai(end)+1);
  end  
    if rext~=1
     rext=rext+DelT_z/(2*rr);
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
    %Delta1=0;
    Delta_gr=anir(istr)/rr^2;
    if igr_app==-1
     Delta_gr=0;
    end
    if abs(Delta_gr)>0 
     if ick==1
      'emmemix', keyboard
     end
    end
    if ifp>0
       disp('[spessore nint next]')
       [dos rint rext]
       disp('[a n]')
       faid=[fai fai(end)+1];
       [[aitr(fai); 0] riv(faid)]
       pausak
    end

%%%%%%%%% a)

        Kost=0;
%        KostV=0;
        Kostz=0;
        if igr_app==1
         ifianitr=0;
        else
         ifianitr=1;
        end
        if icomp_anymat==2  %reticolo curvo vero
         ifianitr=0;
        end
        
            prtr0=prtr(1);
        
        %for ianu=1:length(fai)   %(ciclo su ogni sezione trasversale omogenea)
            ipri=0;
        for ianu=fai   %(ciclo su ogni sezione trasversale omogenea)
            istr_ret_curv=0;
            at0=aitr(ianu);
            bt0=bitr(ianu);
            pt0=pitr(ianu);
%            prtr0=prtr(ianu);
            sh0=abs(shtr(ianu));
            
            if iSET==1
             ' cont', keyboard
            end
            if sh0==7
             sh01=1;
             
% parte nuova per anisotropia specchio curvo-reticolo
              if ifianitr==0
               fimin=find(prtr<0); 
               ivalany=0;
               if length(fimin)==1 
                fiani=find(rad_any==at0);
                 if length(fiani)==1
                  ivalany=1;
                 end 
               end
              end 
              
              if ifianitr==0 & ivalany==1
               ifianitr=1;
               ick=ickset;

               if length(fimin)>0 
               igra=1;
                if icomp_anymat==2 %reticolo curvo vero
                 istr_ret_curv=1;
                 dialef=diag(pes);
                end 
                

%               'anisotropia specchio'
%                 if icomp_anymat==1  %reticolo curvo appr
%' quo ret', keyboard
                if icomp_anymat==2   %reticolo curvo vero

                  if prtr(fimin)==-1
                   Kan_gr=KAN(:,:,fiani)*diag(pes);
                   iset_pri=1;
                  elseif abs(prtr(fimin))>=2
                   iset_pri=0;
                   if length(find(aitr>0))==1
                    Kan_gr=KAN(:,:,1)*dialef-KAN(:,:,fiani)*diag(pes);     % K_infinito-K_limitato
                   else 
                    fianiu=find(rad_any==aitr(ianu+1));
                    Kan_gr=(KAN(:,:,fianiu)-KAN(:,:,fiani))*diag(pes);     % K_infinito-K_limitato
                   end 
                  end 
                  if iztm==1   %reticolo curvo vero
                   if prtr(fimin)==-1
                    Kanz_gr=KANz(:,:,fiani)*diag(pes);
                   elseif abs(prtr(fimin))>=2
                    if length(find(aitr>0))==1
                     Kanz_gr=KANz(:,:,1)*dialef-KANz(:,:,fiani)*diag(pes);     % K_infinito-K_limitato
                    else 
                     fianiu=find(rad_any==aitr(ianu+1));
                     Kanz_gr=(KANz(:,:,fianiu)-KANz(:,:,fiani))*diag(pes);     % K_infinito-K_limitato
                    end                   
                  end
                 end 

                 end
               end %ifianitr
               
               if icomp_anymat==1  %reticolo curvo vero
                if Psalva.orien==0
                 seg_ret=1;
                else
                 seg_ret=-1;
                end
                rintz=nitr(ianu)+seg_ret*ani_gr/2;  % eps_z = eps_x
               end 
                if ick==1              
                  at0
%                  'sh0 anisotropo', keyboard
                end  
                if istr==10700
%                                  'sh0 anisotropo', keyboard
                  'sh0 anisotropo', keyboard

                end
               
              end
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
             depepfz=-rr^2/riv(ianu)^2+rr^2/rext^2;
             if sh0==6 & igraef_new==0
              rext=sqrt(n_me);
%              depepf=(riv(ianu)^2-riv(end)^2)/rr^2;
%              depepfz=-rr^2/riv(ianu)^2+rr^2/rext^2;              
%              'sistemo gli indici', keyboard
             end
%             depepf=riv(ianu)^2/rr^2-1;
%             depepfz=1-rr^2/riv(ianu)^2;
%'passo??', keyboard
               if sh01==11
                if SPP==0
                 rext=riv(1);
                 depepf=(riv(2)^2-riv(1)^2)/rr^2
                 depepfz=-rr^2/riv(2)^2+rr^2/rext^2
                else
                 rext=riv(length(fai)+1);
                 depepfV=(riv(1)^2-rext^2)/rr^2
                 depepfzV=-rr^2/riv(1)^2+rr^2/rext^2                 
                end
%'rext SPP =11', keyboard

               end

            else
             depepf=(riv(ianu)^2-riv(ianu+1)^2)/rr^2;
             depepfz=-rr^2/riv(ianu)^2+rr^2/rext^2;
%		depepf, pausak
%		keyboard
            end
%		depepf, pausak

             ipriv=ipri;
            if sh01==1
             ipri=find(aytot(:,sh0)==at0);
%             if length(ipri)>1
%              ipri=ipri(1);
%             end
            elseif sh0>1 & sh0<4 | sh0==11
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

            if istr_ret_curv==1
             iama=[];
             PUrid=[];
             if iset_pri==1
              ipri=0;
             end 
%             ' ret vero', keyboard
%             ' ret vero', keyboard
             depepf_gr=(ngra(1)^2-ngra(2)^2)/rr^2;
             depepfz_gr=-rr^2/ngra(1)^2+rr^2/ngra(2)^2;

             Kr=depepf_gr*Kan_gr;
             if ick==1
              'Kost reticolo', pausak
             end 
             Kost=Kost+Kr;
             if iztm==1
              Krz=depepfz_gr*Kanz_gr;
              Kostz=Kostz+Krz;
             end
             if ick==1 
%              ' ret vero', keyboard
               icontrolla=1;
             end
            end 
        if icontrolla==1 
        ' [istr prtr aitr real(nitr)]'
         [istr prtr aitr real(nitr')], pausak
        end
            
%           ' selezione e_mix', keyboard
            Kr=0;
%            istr
%            ' wmme mix', pausak
           if sh0==6 & ired_ret==2
            dpes1=Par{1}.dpes1;
            dpes2=Par{1}.dpes2;
            Idelta=Par{2}.Idelta;
            Ideltaz=Par{2}.Ideltaz;
            Madd0=Par{2}.Madd0;
            Madd0z=Par{2}.Maddz0;
            be=Par{2}.be;
           end 
           if sh0==6 & (igraef_new==2 | igraef_new==0)
            ipri=0;
           end
%           icondSPP=
           if ipri~=0 & sh0~=11
             iama=Iacc{sh0};
             PUrid=PUrD{sh0};
%             iama=[];
%             pausak
%'qui KOS loop', keyboard
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
               Kvi(:,:,ipri)=KOS(:,:,ipri);
               if ianu>1 
                if shtr(ianu-1)==11
                 at0p=aitr(ianu-1);
                 ipri0=find(aytot(:,sh0)==at0p);                 
                 Kr=dpes1*(KOS(:,:,ipri)-KOS(:,:,ipri0))*dpes2*depepf;
                end
               end
%              KrV=KOS(:,:,ipri)*depepf;

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
             'lateral size: '
              ianu
             ' ferma: at0 bt0 pt0 sh0 ipri = ',[at0 bt0 pt0 sh0 ipri], keyboard
            end

             if ick==1
              'Kost normale', pausak
             end
            Kost=Kost+Kr;
%            KostV=KostV+KrV;
            


            
            if iztm==1 & ipri~=0
              KOS=Kosz{sh0};
             iama=Iacc{sh0};
             PUrid=PUrD{sh0};
%             iama=[];
%             pausak
              if  length(size(KOS))==3
%               Krz=dpes1*reshape(KOS(:,:,ipri),si2)*dpes2*depepfz;
                if ianu>1
                 Krz=dpes1*(KOS(:,:,ipri)-KOS(:,:,ipri-1))*dpes2*depepfz;
                else
                 Krz=dpes1*KOS(:,:,ipri)*dpes2*depepfz;
                end 
               if ianu>1 
                if shtr(ianu-1)==11
                 at0p=aitr(ianu-1);
                 ipri0=find(aytot(:,sh0)==at0p);                 
                 Krz=dpes1*(KOS(:,:,ipri)-KOS(:,:,ipri0))*dpes2*depepfz;
                end
               end                
              else
%               Krz=KOS*diag(pes)*depepfz;
               Krz=dpes1*KOS*dpes2*depepfz;
              end
              clear KOS

              Kostz=Kostz+Krz;
            end

            if ifp>=1
             'dopo Kost e Kostz '
              ianu
              pausak
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
        if icontrolla==1 
        ' [istr prtr aitr real(nitr)]'
         [istr prtr aitr real(nitr')]
          'ret vero: controllo matrice totale', pausak
        end
        if iztm==1
%           KOt=ck*(Idelta*depepo+Delta*Kan+Delta_gr*Kan_gr+Delta1*Kan1+Kost+KTep);
           KOt=ck*(IDE+ianisotropia*(Delta*Kan+Delta_gr*Kan_gr+Delta1*Kan1)+Kost+KTep);
           % Ks=Delta1*Kan1;
           if igra==1 & ick==1
%            'verifica tutto', pausak
           end
%           if Delta_gr~=0fi
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
           KOt=ck*(IDE+ianisotropia*(Delta_gr*Kan_gr+Delta*Kan+Delta1*Kan1)+Kost+KTep);
           P=[KOt KOt; -KOt -KOt];
        end
       if ifp==1
        [istr; aitr], 'emme_mix', pausak
       end
 if istocont==1
       ' Emme P', keyboard
 end       
if abs(Delta_gr)>0
 iama=[];
% 'Delta', keyboard
end
if sh0==6
% 'Controllo ret vero', keyboard
end
%        Mo=diag([Belo; -Belo])+P;

        Belod=([Belo; -Belo]);
        fidia=find(diag(Idelta)==0);
        Belod(fidia)=0;
% '        passo qui', keyboard
        P=diag(Belod)+P;

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
%'ilayfast', keyboard
          if ilayfast(istr)==0 | (Ps.tuttefreq==1 & sh0==11)

%          if (sh0==6 & igraef_new==2)==1

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
            if iSET==1
             ' SET Oo 1', keyboard
            end 
            
            if sh0==6 & igraef_new==2
%             'esatto ', keyboard
              pol
%              pausak
if length(find(pvet(1)==pol))==1  

tuttefreq=1;
tuttefreq=0;
if isfield(Ps,'tuttefreq')==1
 tuttefreq=Ps.tuttefreq;
end

global iORTA
%iORTA=1
%' dopo iOrta', keyboard

  if tuttefreq==1
    condiz_ifr=1;
    lar=1/kcav*2*pi*rr;
  end  
  if tuttefreq==0
   if ifr==1
    condiz_ifr=1;
    dLa=diff(Dlam_mod(1:2))/2000;  
    lar=1/kcav*2*pi*rr+dLa;
   else
    condiz_ifr=0;
   end
  end   
 if iORTA==0
  condiz_ifr=0;
 end
  if condiz_ifr==1 & iFie==0
%                isao=1;
%               if ifr~=0
%                isao=0;
%               'calcolo ORTA',             keyboard
               'calcolo ORTA',            
            
            
%    save corta
    mbvero=mm+[0:pasnu:numodiacc*pasnu];
               ifpT=ifp;
%                'lar Orta ',  keyboard

   
               fir=find(shavet==6);
               nr1=nv(fir,1);
               nr2=nv(fir,2);               
               r_in=nv(fir-1,1);
               r_out=nv(fir+1,1);

%               'qui Optica ',  keyboard

	if iBEL==1963
           [Oo1,Oo2,G1,G2]=Teq1_modif(KK,lar,dos,period,DC,nr1,nr2,r_in,r_out,rr,mbvero,ifpT,segem);               
	else
% prima di sett 7 2015	
%           [Oo1,Oo2]=Teq1_modif2013(KK,lar,dos,par_grat,rr,mbvero,ifpT,segem,icriti);               
        isca=1;
          if mm==mmvet(1) 
%               'calcolo ORTA',             
           [Oo1,Oo2,Tort]=Teq1_VORTEX(KK,lar,dos,par_grat,rr,mbvero,ifpT,segem,icriti,isca);               
           Torts{ifr}=Tort;
          else
               'calcolo Solo Proiezioni',             
           Tort=Torts{ifr};
           [Oo1,Oo2]=Teq1_VORTEX(KK,lar,dos,par_grat,rr,mbvero,ifpT,segem,icriti,isca,Tort);               
          end
	end
%               'Salvo Orta', keyboard
      
      if mm==1
%     		OSAV1.o1=Oo1(Pusta,Pusta);
%     		OSAV1.o2=Oo2(Pusta,Pusta);
     		OSAV1.o1=Oo1;
     		OSAV1.o2=Oo2;     		
%     		OSAV1.o3=Oo3;
%     		OSAV1.o4=Oo4;

        else
%     		OSAV0.o1=Oo1(Pusta,Pusta);
%     		OSAV0.o2=Oo2(Pusta,Pusta);   
     		OSAV0.o1=Oo1;
     		OSAV0.o2=Oo2;        		
%     		OSAV0.o3=Oo3;
%     		OSAV0.o4=Oo4;        		
     		
        end             

 end       % calcolo orta
end        % pol 

            if mm==1
             if iORTA==0 & ifr==1 & pvet(1)==pol
%              'carico Orta 1', pausak
%               eval(['load ',Dsav,'\saOrta1 ']);     		
%               load saOrta1
             end
             OSAV=OSAV1;
            else
             
             if iORTA==0 & ifr==1 & pvet(1)==pol
%              'carico Orta 0', pausak
             
%               eval(['load ',Dsav,'\saOrta0 ']);     		
%               load saOrta0
             end
             OSAV=OSAV0;
            end
               if pol==1
                Oo=OSAV.o1;
%                Oo=Oo2;
%                Ga1=Ga1*0;
%               ' qui pol 1', pausak
               elseif pol==-1
%                Ga1=Ga1*0;
%                Oo=Oo2;
                Oo=OSAV.o2;
                
%                Oo=Oo1;
%               ' qui pol 2', pausak
               elseif pol==0
                Oo1=OSAV.o1;               
                Oo2=OSAV.o2;
%                Oo3=OSAV.o3;               
%                Oo4=OSAV.o4;                
                pupo=length(Oo1)/2;
                p11=1:pupo;
                p22=p11+pupo;
                Ze=zeros(pupo,pupo);
                M11=[Oo1(p11,p11) Ze ; ...
                       Ze Oo2(p11,p11) ];  
                M12=[Oo1(p11,p22) Ze ; ...
                       Ze Oo2(p11,p22) ];       
                M22=[Oo1(p22,p22) Ze ; ...
                       Ze Oo2(p22,p22) ];  
                M21=[Oo1(p22,p11) Ze ; ...
                       Ze Oo2(p22,p11) ];  

   %             M11=[Oo1(p11,p11) Oo3(p11,p11) ; ...
   %                    Oo4(p11,p11) Oo2(p11,p11) ];  
   %             M12=[Oo1(p11,p22) Oo3(p11,p22)  ; ...
   %                    Oo4(p11,p22)  Oo2(p11,p22) ];       
   %             M22=[Oo1(p22,p22) Oo3(p22,p22)  ; ...
   %                    Oo4(p22,p22)  Oo2(p22,p22) ];  
   %             M21=[Oo1(p22,p11) Oo3(p22,p11) ; ...
   %                    Oo4(p22,p11)  Oo2(p22,p11) ];  
   %                    
                 OoNR=[M11 M12;
                       M21 M22];
                       
               Oo=[M11(Pus0,Pus0) M12(Pus0,Pus0);
                       M21(Pus0,Pus0) M22(Pus0,Pus0)];
%                   Oo=OoNR;
                       Oogr=Oo;
%                       'ver grap TEMP', keyboard
               end       
               if abs(pol)==1
                                pupo=length(Oo)/2;
                                PuT=[Pus0 Pus0+pupo];
                                Oo=Oo(PuT,PuT);
               end
            else
           if isnan(P(1))==1
           %if ifp>0 & ick==1
           'Oo NAN', keyboard
           'Oo', keyboard
           'Oo', keyboard
           'Oo', keyboard
           end
%'GEREP', keyboard
             parmatcrit=(abs(depepo)+abs(depepf))*dos*5;
             if parmatcrit>1
              parmatcrit
              iautm=1
              if ifp==-10
               keyboard
              end 
             end
%           'prima di expM', keyboard      
              if length(find(shtr(fai)==1))==length(fai)
               if ~exist('FiMATo')
                subsemo
               end              
               ifpla=0;
               iama.Fi=FiMATo;
               iautm=0;
%               ifpla=-100;               
              else
               iama=[];
              end           
%             'qui Vor', keyboard
              if length(find(shtr==11))==0
               Oo=expm_mioVor(P,ifpla,iama,Pust,PUrid,iautm);
              else
               if SPP==1
	        OoSPPmet               
               else
                Oo=expm_mioVor(P,ifpla,iama,Pust,PUrid,iautm);
               end
              end 
%             'dopo Vor', keyboard
             if sh0==6
              OogrBW=Oo;
              'calcolato OogrBW'
              keyboard
             end
             
%               if sh0==11 & SPP==1
%%'entro in OoElica', keyboard
% 	        OoSPP
%%	        'find OoElica', keyboard 	        
%  	       end
  	       
            end
           end
           if isnan(Oo(1))==1
           %if ifp>0 & ick==1
           'Oo NAN', keyboard
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
           istri=istrv(istr);
           icousa=icousav(istri);
%           if ifr==2
%            'qui non sem', pausak
%            keyboard
%           end
            if iSET==1
             ' SET Oo 2', keyboard
            end 
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
            if iSET==1
             ' SET Oo 3', keyboard
            end 
            Oo=Tstof+Tstof*P(Pust,Pust)*(freq-fre_camp(pufreq));
           end

          end  %ilayfast





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Se la struttura trasversa e' banale   (omogenea)

 else % aitr  ai
   ifpla=1;
   if ifp>=2
      disp('sono entrata della parte senza struttura trasversa'), pausak
   end
%      disp('sono entrata della parte senza struttura trasversa'), keyboard

%   dos=Li(istr);
   Delta=ani(istr)/rr^2;
   Delta1=ani1(istr)/rr^2;
%   Delta1=0
   Delta_gr=anir(istr)/rr^2;
   
       if igr_app==-1
        Delta_gr=0;
       end
       if sum(abs([Delta_gr Delta Delta1]))~=0
        ifpla=2;
%        ifpla=0;
       end
        if pasnu==1 & ipolar==0
         ifpla=0;
        end       
     if ifp>=2
       disp('[spessore nint larghezza anisotropia]')
       [dos rint ai(istr) Delta]
       pausak
     end
     if rint~=1
      rintv=rint;
      rint=rint+DelT_z/(2*rr);
     end
   if prtr<0
    Kan_gr=KAN(:,:,1)*dialef;     % K_infinito
                if Psalva.orien==0
                 seg_ret=1;
                else
                 seg_ret=-1;
                end
                rintz=nitr(1)+seg_ret*ani_gr/2;  % eps_z = eps_x    
       if ick==1
        'a=0', keyboard
       end 
   end

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
          
%          P=[-(KOt+KOz)  -(KOt-KOz); 
%              (KOt-KOz)      (KOt+KOz)];
%          ' P ', keyboard
       else
          KOt=ck*(Idelta*depe+Delta*Kan+Delta_gr*Kan_gr+Delta1*Kan1+KTep);
          P=[KOt KOt; -KOt -KOt];
%          Pv=[KOt KOt; -KOt -KOt];
       end
%       'Emme piano', keyboard
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
if rint==1 & dos>1
%  'airgap', keyboard
end
            P_sav=P;
            P=diag(Mod)+P;
            
%            P=diag(Mod)+P(Pust);
%'expm_mio 2'
%pausak
%if rint==1 & dos>1
% 'verifica tutto', keyboard
if dos>1
 ilayfastistr=0;
end 
%end
      %      ifpla
      %      if ifpla==2
      %       ifpla=0;
      %      end
          if ilayfastistr==0
%   	 	if ifp==-10
%   			istr, 'Calcolo completo Trasm'       
%   		end
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
            if iSET==1
             ' SET Oo 4', keyboard
            end 

            expm_miu(P,ifpla,iama,Pust);
%            size(Oo)
%            ' emme ver', pausak
           else
            if max(max(abs(KTep)))>0 
             ifpla=0;
             iama.Fi=FiMAT;
             iautm=0;
             ifpla=-100;
            else
              if length(find(shtr(fai)==1))==length(fai)
               if ~exist('FiMATo')
                subsemo
               end              
               ifpla=0;
               iama.Fi=FiMATo;
               iautm=0;
               ifpla=-100;               
              else
               iama=[];
              end
            end
%'qui pla', keyboard            

%tic
            Oo=expm_mioVor(P,ifpla,iama,Pust,PUrid,iautm);

%            'tempo nuovo'
%toc     
%tic
%            Oov=expm_mio(P,0,[],Pust,PUrid);
%            'tempo vecchiio'
%toc     
%map(log10(abs(Oo-Oov+1e-18)))
%             ' SET Oo 5', pausak
if istr==88
 %'88', keyboard
end
            if iSET>=1
%             ' SET Oo 5', pausak
%             load saf
%             load saO
             if pol==1
%              Oo=Ttots;
%              Oo=Oo1;
              Oo1=expm(P);
              ' fine Oo1', keyboard
             else
%              Oo=Ttotc;
%              Oo=Oo2;
              Oo2=expm(P);
              if sh0==6
              r1=nv0(4,1);
              r2=nv0(4,2);              
              DC=t1/(t1+t2);
              end
              'fine: salvare', keyboard
             end
            end            
%           ' SET Oo 5', keyboard

             if rint==1
%               ' SET Oo 5 dopo', keyboard  
             end
            if isnan(Oo(1))==1
           %if ifp>0 & ick==1
           'Oo NAN', keyboard
           'Oo', keyboard
           'Oo', keyboard
           'Oo', keyboard
            end
           end

           freqs=freq;

          else  %ilayfast

           icousa=icousav(istr);
           if iTsav==0
%            Oo=Tstof{icousa}+Tstof{icousa}*Mo(Pust,Pust)*(freq-fre_camp(pufreq));
%            ' SET Oo 6', keyboard
            if iSET==1
             ' SET Oo 6', keyboard
            end 
            Oo=Tstof{icousa}+Tstof{icousa}*P(Pust,Pust)*(freq-fre_camp(pufreq));
           else
            if ispeed==1
             eval([' load ',Dsav,'\',nTstof,num2str(icousa)]);
            end
%            Oo=Tstof+Tstof*Mo(Pust,Pust)*(freq-fre_camp(pufreq));
            if iSET==1
             ' SET Oo 7', keyboard
            end 
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
            if iSET==1
             ' SET Oo 8', keyboard
            end 
    Oo=IdeOon;
    
   end

  end  %ai
  
  if istr==88
  % 'reset Oo'
  % Oo=eye(size(Oo));
  end 
%     istr
%     figure, plot(diag(abs(Oo))), pausak
if istr==3 & dovesono==1 & sh0==3
' dopo ', pausak
end
if ifp>0 | ifp==-11
'istr=', istr
' qui emme_ult',% keyboard
pausak
end
if length(Oo)~=length(Tc)
%'qui pro', keyboard
end
%'qui pro scala', keyboard

     if imem==0
      Tc=Oo*Tc;
     else
%     istr
%     dovesono
%' qui dopo emme_ult', keyboard
      Tc=prodmat(Oo,Tc,ifpla,Pust);
     end
if ifp>0 | ifp==-11
'istr=', istr
' qui dopo emme_ult',% keyboard
pausak
end

%  if ifp>=2 & length(Oo)>5 | (rint==1 & dos>1)
%  if ifp>=2 & length(Oo)>5 | (sh0==6)
  if ifp>=2 & length(Oo)>5 
   figure, surf(abs(Oo)), shading('interp'), view(0,90),
   title(' Oo ')
   colorbar, pausak
  end

else
%            ' SET Oo 9', keyboard
 Oo=IdeOon;

end  %layskip

if ifp~=-4 & ick==1
[istr]
'emme', pausak
end

if ifiez==0
if layskip(istr)==0 & ilaymem(istr)==1 & (ifr==1 | ilin==0)
   if pol==pvet(1)
    icoustor=icoustor+1;
%    icoustor=istr;
    icousav(icoustor)=icoustor;
    istrv(istr)=icoustor;
   else
    istri=istrv(istr);
    icoustor=icousav(istri);
   end
   %if dovesono==2
   %'fine'
   %[icoustor istr]
   %pausak
   %end
   
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
   Prga_ret

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
if dovesono==3 & istr>216
%if ani_gr~=0
istr

' % dopo calcolo M ', keyboard
end
if ifp>0
' % dopo calcolo M '
 if ifp>2
  keyboard
 end
end

if icrit(istr)==1
% icrit(istr)=1;
 Mcrit{istr,icpo}=P;
  if sh0==6 & igraef_new==2
    Mcrit{istr,icpo}=Oo;
  end
% 'Mcrit', pausak
end 

if  istr==88 & igraef_new==1 & lp1==1
 if pol==1
  Oo1=Oo;
 else
  Oo2=Oo;
 end
 if exist('Oo1')*exist('Oo2')==1
  save saBW Oo1 Oo2
 end
%'Any0', keyboard
end
%if istr==5 & dovesono==1
if istr==4 & dovesono==10
%if istr==-7 & dovesono==1
 istr
 dos
% ifp=3;
 'dove', keyboard
end

%if max(max(abs(KTep)))>0
% ' controlla Temp', pausak
%end
if sh0==6
%' fine  emme scala', keyboard
end
