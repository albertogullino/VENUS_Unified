
pufreq=ifr-1;
ilin=0;
istan=0;
iztm=1;


istr00=10+2000;

PerStop=1;
%PerStop=1000;
%PerStop=dovesono-1;

istrSTOP=9*PerStop;


istrSTOP=[7 9 16];
istrSTOP=[9 16];
istrSTOP=[16];
istrSTOP=[13]*PerStop;
istrSTOP=[4]*PerStop;  % anisotropia laterale
%istrSTOP=[7]*PerStop;  % grating
%istrSTOP=[14]*1000;

DSO=1;
if istr==istrSTOP
 'contollo tutto', keyboard
end

%'entro emme', keyboard




ifp=-4;

iStrCont=5;
%iStrCont=1;
iCont=0;
Kr=0;
Krz=0;
DelsuEpst=0;
DelsuEpsz=0;


icontrolla=0;
if min(abs(istr-istrSTOP))==0 & dovesono==1
 ifp=3
 icontrolla=1;
end

ideold=0;

AnyCirc=0;
iDelta=1;
%kcav=1.074717021057574e+001


global Dsav
Oo=IdeOon;
iSET=0;
if istr>=istr00
 iSET=1;
end
if iSET==1
'inizion emme', keyboard
end
sh0=-1000;
aiduM=1;

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
   igraef_new=-1;
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
         elseif ifiez==1
          dos=Liz;
          ilayfastistr=0;
          layskipistr=0;
         elseif ifiez==2
          ilayfastistr=0;
          layskipistr=0;
         end

if istr>=8 & dovesono==10
 'controllo Any', keyboard
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



PUrid=0;

if length(find(shavet==6))>0
 iproga=1;
else
 iproga=0;
end

 iproga=0;

[Tfas,Ffas]=eltime(Tfas,Ffas,istr);

if ifp>=2
%disp(' entro emmestru: istr ')
%istr
end

%keyboard

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
    'salvato'
    keyboard

   end
   ifplatot=1*isem;
   T=IdeOon;
  end
end



%if layskip(istr)==0
if layskipistr==0


 if dovesono~=3
  prtr=prosi(:,istr);
  aitr=ai(:,istr);
  pitr=pai(:,istr);
  fier=find(aitr>0);
  aierr=aitr(fier);
  if length(find(diff(aierr)<0))>0 & dos>0
   ' in Emme_any0_11'
   ' File errato !!!!  raggi in ordine errato ', keyboard
  end
  bitr=bi(:,istr);
  if aitr>0
   pitr=pai(:,istr);
  end
  shtr=shai(:,istr);
  


  if ity==1
   tytr=tyari(:,istr);
  end

  nitr=ni(:,istr);
 else
 % regione attiva
  nitr=ni;
  aitr=ai;
  bitr=bi;
  pitr=pai;
  shtr=shai;
  prtr=1;
 
 end
  
  fi4=find(shtr==4);
  if length(fi4)==0      % anisotr. piana
   fi4=1;
  end
  if length(fi4)>0
   if length(fi4)>1
    findEQ=find(nitr==nitr(1));
    fi4=findEQ;
   end 
%   'fi4', keyboard
    if isfield(Ps,'nOrd')
     nOrd=Ps.nOrd;
     nExt=Ps.nExt;
     if pol==1
      fiEx=find(nitr(2:end)==nExt);
      nitr(fiEx+1)=nOrd;
     end
    end 
  end  
  
  if length(find(imag(nitr)<-1))
%   ifp=3
%   keyboard
  end
  aniloc=ani(istr);
  Aniloc=0;
  Anilocz=0;
  naverage=nitr(1);
  clear nz
if abs(imag(aniloc))>0 & iany==3
  
%    'ani', keyboard
if ~exist('iany_CIRC')
 iany_CIRC=0;
end
  if iany_CIRC==0  
   nord=imag(aniloc);
   nextrM=-real(aniloc);
   if nextrM<naverage
    if abs(nextrM-naverage)<.1
     fi=find(nitr==naverage);
     naverage=nextrM;
     nitr(fi)=naverage;
    else

    'errore cristallo liquido: aggiusta file str', keyboard
    end
   end
   if nord*nextrM<0
    'errore cristallo liquido: aggiusta file str', keyboard
   end   
%   nextrM=Par.nextr;
   DelLC=nextrM-nord;
   
   for ianytr=fi4'
   nextr=nitr(ianytr);
   ndelta=0;
   if DelLC~=0
   teta=acos((nextr-nord)/DelLC);
   nz=nord+DelLC*sin(teta);
   ey=nz^2;
%   ey=nextr^2;
   naverage=sqrt((nextr.^2+nord.^2)/2);
   ndelta=(nord.^2-nextr.^2);
   Aniloc(ianytr)=ndelta/rr^2;
   %Aniloc=1;
   Anilocz(ianytr)=1-rr^2/ey;
   end
 
    nitr(ianytr)=naverage;
   end %
   
   if length(find(shtr==-8))==0
    ani(istr)=ndelta;
   else
    ani(istr)=0;      % anisotropia confinata � gi� in sh=-8
   end
%   'qui ani', keyboard
   ifpla=0;
   if ifr==1 & ifp==-10
    'ANY ver', keyboard
   end
  else %CIRC
   iDelta=0;
   nextr=nitr(1);
   nitr(1)=rr;
   nord=aniloc;
   nextrM=Par.nextr;
   DelLC=nextrM-nord;
   teta=acos((nextr-nord)/DelLC);
   nz=nord+DelLC*sin(teta);
   ey=nz^2;
   n_ro=nord;
   n_fi=nextr;
   DelTE=2*(n_fi.^2/rr^2-1);
   DelTM=2*(n_ro.^2/rr^2-1);
   IdeA=eye(length(KKv));
   IdeAe=IdeA;
   IdeAm=IdeA;
   if mm==0
    if pol==1
     IdeAe(1:lk,1:lk)=0;
    elseif pol==-1
     IdeAm(1:lk,1:lk)=0;
    else
    'errore pol anis. CIRC', keyboard
    end
   end
   
   DeMat=[IdeAe*DelTE IdeA*0; IdeA*0 DelTM*IdeAm];
   AnyCirc=DeMat.*IdeltaaC;

   'ANY Circ', keyboard 
  
  end

 end
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


%    disp('  istr, ifpla, ifplaco, ifplatot ')
%   [istr, ifpla, ifplaco, ifplatot],
%   pausak

%'qui fai', keyboard
  fai=find(aitr~=0)';

  if length(fai)==0
   fai=1;
  end 
  

  fi8=find(shtr==-8);
  if length(fi8)>0
%   'QIO', keyboard
   if length(nitr)==2
    nitr=[nitr; 0];
   end
%   riv(fi8)=rr;

%   'fait', keyboard
%   fai=1;
  end
  rint=nitr(fai(1));
  rintz=rint;

% parte per profilo di perdite trasversale


      pim=imag(rint^2);

 DelT_z=0;
 if igau==4 
  KmeT=mean(mean(KTe));
 else
  KmeT=0;
 end 
 if KmeT~=0

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
    riv=nitr; 
    rivz=nitr; 
    if exist('nz')
     rivz(fi4')=nz;
    end 
    
  if length(fi8)>0
%   riv(fi8)=rr;
     if iBWnew==0 
               RA=radii.array{istr+1};
               pgra=RA{13};
               igra=pgra.igra;

%              ' da verificare! ', keyboard
		   g_n1=par_grat.n1;
		   g_n2=par_grat.n2;
		   NModes=par_grat.NModi;
                   Period=par_grat.period*1000;
                   DC=par_grat.DC;
		   lambda1000=lambda*1000;
		   thickness=Li(istr)*1000;
		   if istr<fiQW
		    nin=ni(1,istr-1);
		    nout=ni(1,istr+1);
		   else
		    nin=ni(1,istr+1);
		    nout=ni(1,istr-1);		   
		   end                   
		   
		   [neff,nBW,flagt] = f_EffectiveIndex(Period,lambda1000,DC,nin,nout,g_n1,g_n2,thickness,NModes);
		   if igraef_new==0
		    n_pa=nBW(1);
		    n_ve=nBW(2);
		   else
		    n_pa=neff(1);
		    n_ve=neff(2);  
		   end

	           if pol==1
	            ngr=n_ve;
	           else
	            ngr=n_pa;
	           end
	           figr=find(shtr==-8);
	           riv(figr)=ngr;
	
%	           'qui riv', keyboard

    end % iBWnew

   
   
   
  end    
  if length(fi8)>0
   sh0=8;
%   ' riv vedo ', keyboard
  end 
%   istr
%   fai'
%   pausak

  end
%  if istr==8
%  ' istr=32', keyboard
%  end
  if ibrea==1
   ' qui' , keyboard
  end
%  if aidu==0
  if aidu<0 & sh0~=8
   dos=0;
  end
  if aidu<0 & sh0==8
   aiduM=aidu;
   aidu=abs(aidu);
  end  
%%%%%db%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Se lo strato presenta una struttura trasverale non banale
% if aidu~=0 & Li(istr)~=0

% if aidu~=0 & dos~=0
%istr
%keyboard


 if aidu>0 & dos~=0
 
   if length(riv)==1
     rext=riv(end);
   else
     rext=riv(fai(end)+1);
%     rext=riv(length(fai)+1);
   end  
   
     if rext~=1
      rext=rext+DelT_z/(2*rr);
    end
 
     if ifp>=2
        istr
        display('sono entrata della parte con struttura trasversa'), pausak
     end

%    dos=Li(istr);

%    Delta=0;  % anistropia confinata trattata a parte da Aniloc
%    Delta1=0;  % anistropia confinata trattata a parte da Aniloc
%   Delta_gr=0;  % anistropia confinata trattata a parte da Aniloc
% if ani(istr)==0
%    Delta_gr=anir(istr)/rr^2;
% end    
   if length(Aniloc)==1 
    iDelta=0;
   end
   
%  'qui Delta', keyboard 
   Delta=iDelta*ani(istr)/rr^2;
   Delta1=ani1(istr)/rr^2;
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
        Kost=K_antiguiding;
        Kostz=K_antiguidingZ;
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
         
         Delz=(1-rr^2./rivz.^2);
         Delt=(riv.^2/rr^2-1);
	     for ktr=1:fai(end)
	      DelsuEpsz(ktr)=Delz(ktr)-Delz(ktr+1);
	      DelsuEpst(ktr)=Delt(ktr)-Delt(ktr+1);
	     end
	    

	   if length(Aniloc)>1 
	    AniSa=Aniloc; 
	    AnizSa=Anilocz; 	   
	    clear Aniloc Anilocz
	     for ktr=fai
	      Aniloc(ktr)=AniSa(ktr)-AniSa(ktr+1);
	      Anilocz(ktr)=AnizSa(ktr)-AnizSa(ktr+1);
	     end	    
	   end
                       if length(fai)>1
%			       'qui Delz', keyboard
		       end         
  
         %-diff(Delz);
         
            ipri=0;
  for ianu=fai   %(ciclo su ogni sezione trasversale omogenea)
%                               if length(find(imag(nitr)<-1))>0
%                                 ianu
%				       'qui Delz', keyboard
%			       end
            Kany3=0;  % anisotropia finita
            ipriANY=[];

     
            
            
            
            
            
            istr_ret_curv=0;
            at0=aitr(ianu);
            bt0=bitr(ianu);
            pt0=pitr(ianu);
%            prtr0=prtr(ianu);
            sh0=abs(shtr(ianu));
            if sh0==8
  %          'sh0', keyboard
            end
            
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
                  if icomp_anymat==2 & iztm==1   %reticolo curvo vero
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
      
             if ideold==1
       	      depepf=(riv(ianu)^2-rext^2)/rr^2;             
              depepfz=-rr^2/riv(ianu)^2+rr^2/rext^2;
             else 
              depepf=DelsuEpst(ianu);
              depepfz=DelsuEpsz(ianu);
             end 
             
%             depepf=riv(ianu)^2/rr^2-1;
%             depepfz=1-rr^2/riv(ianu)^2;
            else
             if ideold==1
       	      depepf=(riv(ianu)^2-rext^2)/rr^2;             
              depepfz=-rr^2/riv(ianu)^2+rr^2/rext^2;
             else 
              depepf=DelsuEpst(ianu);
              depepfz=DelsuEpsz(ianu);
             end             
            end
            

             ipriv=ipri;
            if sh01==1
             ipri=find(aytot(:,sh0)==at0);
%             if length(ipri)>1
%              ipri=ipri(1);
%             end
            elseif sh0>1 & sh0<4
             ipri=find(aytot(:,sh0)==at0 & axtot(:,sh0)==bt0);
            elseif (sh0==4 | sh0==6)
             if sh0==4
     %        'verdDD', keyboard
              ipri=find(aytot(:,sh0)==at0 & axtot(:,sh0)==bt0 & pdtot(:,sh0)==pt0);
             elseif sh0==6 & igraef_new<2
              ipri=find(aytot(:,sh0)==at0 & axtot(:,sh0)==bt0 & pdtot(:,sh0)==pt0);
             end 
%             if length(ipri)>1
%              ipri=ipri(2);
%             end
            elseif sh0==5
             ipri=find(aytot(:,sh0)==at0 & tytot==ty0);
            elseif sh0==8
%            'qui *88', keyboard
%             ipri=1;
             ipri=find(aytot(:,sh0)==at0);             
            end
            if sh0>1 & sh0<6 & iany==3
%             'anisotropia finita', keyboard
              if sh0==4
               at0ANY=1i*at0;
               ipriANY=find(aytot(:,sh0)==at0ANY);
              end % sh0=4         
             
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
        ' FINE ianu: [istr prtr aitr real(nitr)]'
         [istr prtr aitr' real(nitr')], pausak
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
           if sh0==6 & igraef_new==2
            ipri=0;
           end
     
           if ipri~=0
             iama=Iacc{sh0};
             PUrid=PUrD{sh0};
%             iama=[];
%             pausak
             KOS=Kos{sh0};
             KOSz=Kosz{sh0};   
%             'sha=1',             keyboard
          if sh0==8 & aiduM==1
             'sha=8',             keyboard          
               ndo=size(KOS,3);
               [fi8,co8]=find(shavet==-8);
               RA=radii.array{fi8};
               pgra=RA{13};
               igra=pgra.igra;
               
            if igra==1
                KOSd=KOS(:,:,1);
                KOSzd=KOSz(:,:,1);
                fi88=find(shtr==-8);
                aany=aitr(fi88);
                fi888=find(aytot(:,1)==aany);
                KOSc=Kos{1}(:,:,fi888);
                KOSzc=Kosz{1}(:,:,fi888);                
% parte nuova Apr 18 per anisotropia con pi� sezioni trasversali

                if ndo>1 
                 if iBWnew==2
                  KOSd=(KOS(:,:,2)-KOS(:,:,1));
                  KOSzd=(KOSz(:,:,2)-KOSz(:,:,1));
                 else
                  KOSd=(KOS(:,:,2));
                  KOSzd=(KOSz(:,:,2));                 
                 end
                end
                
%                    Kr=dpes1*KOSd*dpes2;
%                    Krz=dpes1*KOSzd*dpes2;   
		if istr==istr00
%              'Any', keyboard
	        end
%              'Any', keyboard	        
		     KOSt=KOSd*Aniloc;
		     KOSz=KOSzd*Anilocz;
		     for ktr=fi88
	              depepf=DelsuEpst(ktr);
	              depepfz=DelsuEpsz(ktr);	     
	              KOSt=KOSt+depepf*KOSc;
	              KOSz=KOSz+depepfz*KOSzc;
	             end 
	             
                    Kr=dpes1*KOSt*dpes2;
                    Krz=dpes1*KOSz*dpes2;

              if iBWnew==0
	
%	           'qui riv', keyboard

		     KOSt=0;
		     KOSz=0;
		     for ktr=figr
	              depepf=DelsuEpst(ktr);
	              depepfz=DelsuEpsz(ktr);	     
	              KOSt=KOSt+depepf*KOSd;
	              KOSz=KOSz+depepfz*KOSzd;
	             end 
                    Kr=dpes1*KOSt*dpes2;
                    Krz=dpes1*KOSzd*dpes2;
                    
%                    'end BW', keyboard
		end % iBWnew
                if exist('Sig_eps')
		    rint=naverage;
                end
             else % igra
% vero Doe ad anelli        NUOVO       
               KOS=sum(KOS(:,:,1:2:end),3)-sum(KOS(:,:,2:2:end),3);   
               KOSz=sum(KOSz(:,:,1:2:end),3)-sum(KOSz(:,:,2:2:end),3);   
               if isfield(pgra,'rDOE') 
                rivi(1)=pgra.n1;
                rivi(2)=pgra.n2;
                Delzi=(1-rr^2./rivi.^2);
                Delti=(rivi.^2/rr^2-1);
   		ktri=1;
		DelsuEpszDOE=Delzi(ktri)-Delzi(ktri+1);
		DelsuEpstDOE=Delti(ktri)-Delti(ktri+1);
	       else
   		ktri=1;	       
		DelsuEpszDOE=DelsuEpsz(ktri);
		DelsuEpstDOE=DelsuEpst(ktri);		
	       end 
        
	       
	              depepf=DelsuEpstDOE;
	              depepfz=DelsuEpszDOE;	     
	              KOSt=depepf*KOS;
	              KOSz=depepfz*KOSz;
	             
                      Kr=dpes1*KOSt*dpes2;
                      Krz=dpes1*KOSz*dpes2;	       
%	                  'qui doe;', keyboard       
	       
             end

            else %sha=8

              depepf=DelsuEpst(ianu);
              depepfz=DelsuEpsz(ianu);

              if sh0==8 
%               'sono qui;', keyboard
               if  length(size(KOS))==3
                Krz=dpes1*KOSz(:,:,ipri)*dpes2*depepfz;
               else
                Krz=dpes1*KOSz*dpes2*depepf;
               end     

              end

              if  length(size(KOS))==3
%              Kr=dpes1*reshape(KOS(:,:,ipri),si2)*dpes2*depepf;
               Kr=dpes1*KOS(:,:,ipri)*dpes2*depepf;
              else
%              Kr=KOS*diag(pes)*depepf;
               Kr=dpes1*KOS*dpes2*depepf;
              end             
              
              if length(ipriANY)==1  % con Anysotr. trasversale finita
              
%              'qio final;ame ANY3', keyboard
                Kr=Kr+dpes1*KOS(:,:,ipriANY)*dpes2*Aniloc(ianu);
              end
              
            end     %sh0=8
            
%' qui gau', keyboard
              if igau==5 
               KTep=dpes1*KTemp*dpes2/rr^2;
    	       if iztm==1
	        KTepz=dpes1*KTempz*dpes2/rr^2;
	       else
	        KTepz=0;
               end
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

     %        if iSET==1
%              'Kost ', pausak
       %      end
            Kost=Kost+Kr;

            if iztm==1 & ipri~=0
               iama=Iacc{sh0};
               PUrid=PUrD{sh0};
              if sh0~=8
               KOSz=Kosz{sh0};
               if  length(size(KOS))==3
                Krz=dpes1*KOSz(:,:,ipri)*dpes2*depepfz;
               else
                Krz=dpes1*KOSz*dpes2*depepfz;
               end   
% non le aggiungo perch� sono ZERO le matrici
%               if length(ipriANY)==1  % con Anysotr. trasversale finita
%              
%              'qio final;ame ANY3', keyboard
%                Krz=Krz+dpes1*KOSz(:,:,ipriANY)*dpes2*Anilocz(ianu);
%               end                 
              end
            

%             iama=[];
%             pausak



              Kostz=Kostz+Krz;
            end

            if ifp>=1
             'dopo Kost e Kostz '
              ianu
              pausak
            end
           if istr==iStrCont
           if iCont==1
 		'Iacc TRAS', keyboard
           end
           end      
              clear KOS            
%                  ' fine ianu', keyboard
        end  %(fine ciclo sulle sezioni trasversali)
        
            if ifp>=2 & prod(size(Kost))>100
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
        if istr==istr00
        'IDE=', keyboard
        end
        ck=-j*kcav*dos;
        Belo=-j*be*kcav*dos;
        if icontrolla==1 
        ' [istr prtr aitr real(nitr)]'
         [istr prtr aitr' real(nitr')]
          'PRIMA di Oo: controllo matrice totale', pausak
        end
        
        
        if iztm==1
%           KOt=ck*(Idelta*depepo+Delta*Kan+Delta_gr*Kan_gr+Delta1*Kan1+Kost+KTep);
           KOt=ck*(IDE+Delta*Kan+Delta_gr*Kan_gr+Delta1*Kan1+Kost+KTep);
           if iSET==1
            'KOT', keyboard
           end
                 if istr==iStrCont
                  if iCont==1
		   'Verifica tutto', keyboard                   
		  end             
		 end             
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
                       if length(fai)>=1
%			       'qui Delz', keyboard
		       end              
%                       if length(find(imag(nitr)<-1))>0
%	   			' qui Kz', keyboard
%		       end

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

       if ifp>=2
         display('matrice K'), pausak
         figure,  surf(log10(abs(KOt)+1e-5)), shading('interp'), view(0,90),
              title(' KOt ')
         colorbar, pausak
       end
       
 if dovesono==1 & istr==1
%  'verifica T', keyboard
 end       
if abs(Delta_gr)>0
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

%          if ilayfast(istr)==0
         if sh0==6 & igraef_new==2
          ilayfastistr=0;
         end 

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

%           end
           if inuo_bas==0
%                       ' verifica miu ', keyboard
%            size(Oo)
            expm_miu(P,ifpla,iama,Pust);
%            size(Oo)
%            ' emme ver', pausak
           else
            if iSET==1
             ' SET Oo 1', keyboard
            end 
            
            if sh0==6 & igraef_new==2
             if ifp==-10
              'esatto ', keyboard
             end 
%              pausak


		if length(find(pvet(1)==pol))==1  

			tuttefreq=1;
			tuttefreq=0;
		 if isfield(Ps,'tuttefreq')==1
		  tuttefreq=Ps.tuttefreq;
		 end

		global iORTA
%iORTA=1
		  if tuttefreq==1
		    condiz_ifr=1;
		    lar=1/kcav*2*pi*rr;
		  end  
		  if tuttefreq==0
		   if ifr==Pufreq0
		    condiz_ifr=1;
		    dLa=diff(Dlam_mod(1:2))/2000;  
		    lar=1/kcav*2*pi*rr-dLa;
		   else
		    condiz_ifr=0;
		   end
		  end   
		 if iORTA==0
		  condiz_ifr=0;
		 end
		  if condiz_ifr==1
%                isao=1;
%               if ifr~=0
%                isao=0;
		   if ifp==-10
	             'calcolo ORTA',             keyboard
		   end          
               'calcolo ORTA',           % keyboard
            
%    save corta
		    mbvero=mm+[0:pasnu:numodiacc*pasnu];
	               ifpT=ifp;
%                'lar Orta ',  keyboard
	               par_gratDE=radii.array{istr+1};
	              if  length(par_gratDE)==1
	                par_grat=par_gratDE{11};
	              else
	   
	               fir=find(shavet==6);
	               
	               nr1=nv(fir,1);
	               nr2=nv(fir,2);               
	               r_in=nv(fir-1,1);
	               r_out=nv(fir+1,1);
	               par_grat.n1=nr1;
	               par_grat.r1=nr1;
	               par_grat.n2=nr2;
	               par_grat.r2=nr2;               
	               par_grat.r_in=r_in;
	               par_grat.r_out=r_out;               
	              end 
	%               'lar ',  keyboard
%    if segem==1
%     ' segem con reticoli deve essere =1, per coerenza segni TETM in Orta (che ha invertito i segni)'
%     keyboard
%    end

			icriti=1;

%                'Prima Orta', keyboard
		           [Oo1,Oo2]=Teq1_modif2016S(KK,lar,dos,par_grat,rr,mbvero,ifpT,segem,icriti);          
%                'Salvo Orta', keyboard           
	           if ifp==-10
	                'Salvo Orta', keyboard
	            end    
      
		      if mm==1
		     		OSAV1.o1=Oo1;
		     		OSAV1.o2=Oo2;
     %		save saOrta1 OSAV1
%                'mm=1', keyboard
                             if iTsav==1
		                eval(['save ',Dsav,'\saOrta1 OSAV1 ']);  
		             end   

		        else
		     		OSAV0.o1=Oo1;
		     		OSAV0.o2=Oo2;           
%     		save saOrta0 OSAV0
%                'mm=0', keyboard
                             if iTsav==1
		                eval(['save ',Dsav,'\saOrta0 OSAV0 ']);     		
                              end
		        end             

	 end       % calcolo orta
end        % pol 

            if mm==1
             if iORTA==0 & ifr==Pufreq0 & pvet(1)==pol
%              'carico Orta 1', pausak
               if iTsav==1
                eval(['load ',Dsav,'\saOrta1 ']);     
               end 
%               load saOrta1
             end
             OSAV=OSAV1;
            else
             
             if iORTA==0 & ifr==Pufreq0 & pvet(1)==pol
%              'carico Orta 0', pausak
              if iTsav==1
               eval(['load ',Dsav,'\saOrta0 ']); 
              end 
%               load saOrta0
             end
             OSAV=OSAV0;
            end
            if pol==1
                Oo=OSAV.o1;
%                Oo=Oo2;
%                Ga1=Ga1*0;
%               ' qui pol 1', pausak
               else
%                Ga1=Ga1*0;
%                Oo=Oo2;
                Oo=OSAV.o2;
                
%                Oo=Oo1;
%               ' qui pol 2', pausak
               end               
            else
           if isnan(P(1))==1
           %if ifp>0 & ick==1
           'Oo NAN', keyboard
           'Oo', keyboard
           'Oo', keyboard
           'Oo', keyboard
           end
%                       ' verifica mio ', keyboard
             Oo=expm_mio(P,ifpla,iama,Pust,PUrid);
%             if isalvaOo==1
%             ' dopo Oo ', keyboard
%             end
            if iSET==1
             ' dopo Oo ', keyboard
            end 
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
           if exist('istopemme')
           if istopemme==1
            'qui non sem', pausak
           end 
           end
%        Oos=Oo; save sa Oos
%        load sa, mapab(Oos-Oo)
%           end

          else  %ilayfast
%           'qui fast', pausak
           icousa=icousav(istr);
           if ifr==2 & istr==1
%            'qui non sem', pausak
           end
            if iSET==1
             ' SET Oo 2', keyboard
            end 
           if iTsav==0
%            Oo=Tstof{icousa}+Tstof{icousa}*Mo(Pust,Pust)*(freq-fre_camp(pufreq));
%            ' SET Oo 2', keyboard
 	    dfreq=(freq-fre_camp(pufreq));
 	    ifpla;
 	    dovesono;
 	    %pausak
 	    if ilin==1
             Oo=Tstof0{icousa}*(IdeOo+P(Pust,Pust)*dfreq);
            else
             Oo=Tstof{icousa}*(IdeOo+P(Pust,Pust)*dfreq);
            end
            
            if ifr>2
% 	     pausak           
 	    end
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
   Delta=iDelta*ani(istr)/rr^2;
   Delta1=ani1(istr)/rr^2;
   Delta_gr=anir(istr)/rr^2;
 %  Delta_gr=-anir(istr)/rr^2;
   SuD=Delta+Delta1+Delta_gr;
   if SuD~=0
    ifpla=0;
    rintz=riv(fai);
%    rintz=sqrt(ey);
   end
       if igr_app==-1
        Delta_gr=0;
       end
       if Delta_gr~=0
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
   if prtr<0 & igraef_new>0
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
       IDE=depe*Idelta;
        if istr==istr00
        'IDE=', keyboard
        end       
       depez=(1-rr^2/rintz^2);
%       ck=-j*kcav_at*dos;
       ck=-j*kcav*dos;
       Beli=ck*be;
       Mod=([Beli; -Beli]);
       if istr==istrSTOP
        'very ANY3', keyboard
       end
       if iztm==1
          KOt=ck*(IDE+AnyCirc+Delta*Kan+Delta_gr*Kan_gr+Delta1*Kan1+KTep);
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
%          ' P ', keyboard
       else
          KOt=ck*(Idelta*depe+AnyCirc+Delta*Kan+Delta_gr*Kan_gr+Delta1*Kan1+KTep);
          P=[KOt KOt; -KOt -KOt];
%          Pv=[KOt KOt; -KOt -KOt];
       end
%       'dopo P', keyboard
          Pk_sav=P;
    if abs(Delta_gr)>0
     if ifr==Pufreq0 & ifp~=-4 & ick==1
      pr=Idelta*depe+Delta_gr*Kan_gr;
     'emmemix 1', keyboard
     end
    end

       if ifp>=2
         display('matrice K'), pausak
         figure,  surf(log10(abs(KOt)+1e-5)), shading('interp'), view(0,90),
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
   	 	if ifp==-10
   			istr, 'Calcolo completo Trasm'       
   		end
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
%            ' verifica miu 2 ', keyboard
            expm_miu(P,ifpla,iama,Pust);
%            size(Oo)
%            ' emme ver', pausak
           else
            if max(max(abs(KTep)))>0 
             ifpla=0;
            end
if length(find(isnan(P)))>0
 'isnan', keyboard
end
%            ' verifica mio 2 ', keyboard
            Oo=expm_mio(P,ifpla,iama,Pust,PUrid);
            
%             ' SET Oo 5', pausak
if istr==istrSTOP
  'Oo vedo', keyboard
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
            ' SET Oo 6', keyboard
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
  
  if istr==4
%   'reset Oo', keyboard
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


if istr==4
 %'CONT OO', keyboard
end

if ifp~=-4 & ick==1
[istr]
'emme', pausak
end

if ifiez==0
%if layskip(istr)==0 & ilaymem(istr)==1 & (ifr==1 | ilin==0)
if layskip(istr)==0 & ilaymem(istr)==1
   if pol==pvet(1)
    icoustor=icoustor+1;
%    icousav(istr)=icoustor;
    icousav(icoustor)=icoustor;
   else
    icoustor=icousav(istr);
   end
   %if dovesono==2
   %'fine'
   %[icoustor istr]
   %pausak
   %end
%' qui', keyboard
    if icoustor~=istr
     'icoustor e istr diversi! : togliere loops da file str'
     [icoustor istr]
%     pausak
    end 
   
   if iTsav==0
    if length(Tc)==1
     Tstor(:,:,icoustor)=diag(ones(size(Pust)));
    else
      if ncop==0
       Tstor(:,:,icoustor)=Tc;
      else
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
% istr,  'Mcrit', pausak
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
if istr<-8 & dovesono==2
%if istr==-7 & dovesono==1
 istr
 dos
% ifp=3;
 'dove', keyboard
end

if max(max(abs(KTep)))>0
% ' controlla Temp', pausak
end

if dovesono==2
%'emme', keyboard
end

if pol==-1 & dovesono==DSO
%'fine emme', keyboard
end

%'fine emme', keyboard