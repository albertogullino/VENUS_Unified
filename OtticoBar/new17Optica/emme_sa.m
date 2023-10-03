%%%%%%%%%%%/%%%%% All'interno di un ciclo su ciascuno strato istr,
% calcolo M
%tic
%istr
iproga=0;


[Tfas,Ffas]=eltime(Tfas,Ffas,istr);

if ifp>=1
disp(' entro emmestru: istr ')
istr
end

if istr>1
  if layskip(istr-1)==0 & ilaymem(istr-1)==0 & ilaymem(istr)==1
   if pol==pvet(1)
    icoustor=icoustor+1;
    icousav(istr-1)=icoustor;
   else
    icoustor=icousav(istr-1);
   end
%   icoustor
%   pausak

   if length(T)==1
    Tstor(:,:,icoustor)=IdeOo;
   else
    Tstor(:,:,icoustor)=T;
   end
%    Tstor(:,:,icoustor)=T;
%    [icoustor istr]
%    disp(' emme_mix iniz: Tstor'), keyboard
   ifplatot=1*isem;
   T=1;
  end
end



%if layskip(istr)==0
if layskip(istr)==0

 if ilayfast(istr)==0

%  iama=iaccvet{istr};
  aitr=ai(:,istr);
  bitr=bi(:,istr);
  pitr=pai(:,istr);
  shtr=shai(:,istr);
  tytr=tyari(:,istr);

  nitr=ni(:,istr);
  dos=Li(istr);
  aniloc=ani(istr);
  aniloc1=ani1(istr);

    if aitr==0 & igau~=4
     if imem==1
      ifpla=1*isem;
      if (aniloc~=0 & iany==1 & ianys<=1) | (aniloc1~=0 & ianys==1 & iany<=1)
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


 if igau==4
  KTep=KTe*diag(pes)/rr^2;
  KTepz=KTez*diag(pes)/rr^2;
 else
  KTep=0;
  KTepz=0;
 end


  rint=nitr(1);
  fai=find(aitr~=0);
  if length(fai)==0
   aidu=0;
   rext=0;
  else
   aidu=aitr(fai(1));
   ini=find(abs(nitr)==0);
   if length(ini)>0
    riv=nitr(1:ini(1)-1);
   else
    riv=nitr;
   end
   rext=riv(end);
  end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Se lo strato presenta una struttura trasverale non banale
 if aidu~=0 & Li(istr)~=0
     if ifp>=1
        display('sono entrata della parte con struttura trasversa'), pausak
     end

    dos=Li(istr);
    Delta=ani(istr)/rr^2;
    Delta1=ani1(istr)/rr^2;
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
            sh0=shtr(ianu);
            ty0=tytr(ianu);
            if sh0>=5
             depepf=(riv(ianu)^2-riv(end)^2)/rr^2;
             depepfz=-rr^2/riv(ianu)^2+rr^2/riv(end)^2;
%             depepf=riv(ianu)^2/rr^2-1;
%             depepfz=1-rr^2/riv(ianu)^2;
            else
             depepf=(riv(ianu)^2-riv(ianu+1)^2)/rr^2;
             depepfz=-rr^2/riv(ianu)^2+rr^2/riv(ianu+1)^2;
            end

            ipri=0;
            if sh0==1
             ipri=find(aytot(:,sh0)==at0);
            elseif sh0>1 & sh0<4
             ipri=find(aytot(:,sh0)==at0 & axtot(:,sh0)==bt0);
            elseif (sh0==4 | sh0==6)
             ipri=find(aytot(:,sh0)==at0 & axtot(:,sh0)==bt0 & pdtot(:,sh0)==pt0);
            elseif sh0==5
             ipri=find(aytot(:,sh0)==at0 & tytot==ty0);
            end

            Kr=0;

            if ipri~=0
%             if  length(size(Kos))==4
%              Kr=reshape(Kos(:,:,ipri,sh0),si2)*diag(pes)*depepf;
%             elseif  length(size(Kos))==3
%              Kr=reshape(Kos(:,:,ipri),si2)*diag(pes)*depepf;
%             else
%              Kr=Kos*diag(pes)*depepf;
%             end
             iama=Iacc{sh0};
             iama=[];
%             pausak
             KOS=Kos{sh0};
%            disp(' emme_mix')
%            keyboard
             if  length(size(KOS))==3
              Kr=reshape(KOS(:,:,ipri),si2)*diag(pes)*depepf;
             else
              Kr=KOS*diag(pes)*depepf;
             end
             clear KOS
            else
             Kr=0;
            end

            Kost=Kost+Kr;
            if iztm==1
%              if  length(size(Kos))==4
%               Krz=reshape(Kosz(:,:,ipri,sh0),si2)*diag(pes)*depepf;
%              elseif  length(size(Kos))==3
%               Krz=reshape(Kosz(:,:,ipri),si2)*diag(pes)*depepf;
%              else
%               Krz=Kosz*diag(pes)*depepf;
%              end
              KOS=Kosz{sh0};
             iama=Iacc{sh0};
             iama=[];
%             pausak
              if  length(size(KOS))==3
               Krz=reshape(KOS(:,:,ipri),si2)*diag(pes)*depepfz;
              else
               Krz=KOS*diag(pes)*depepfz;
              end
              clear KOS

              Kostz=Kostz+Krz;
            end
            if ifp>0
              disp('[ianu depepf]')
              [ianu depepf]
              figure, surf(abs(Kost)), shading('interp'), view(0,90), colorbar, pausak
              disp(' emmestru ox circ')
              pausak
            end
        end  %(fine ciclo sulle sezioni trasversali)

        if ifp>0
           figure,
           surf(abs(Kost)), shading('interp'), view(0,90), colorbar, pausak
        end

        depepo=(rext^2-rr^2)/rr^2;
        ck=-j*kcav*dos;
        Belo=-j*be*kcav*dos;
        if iztm==1
           depepoz=1-rr^2/rext^2;
           KOt=ck*(Idelta*depepo+Delta*Kan+Delta1*Kan1+Kost+KTep);
           KOz=ck*(Ideltaz*depepoz+Kostz+KTepz);

           Ksu=KOt+KOz;
           Kdi=KOt-KOz;
           P=[Ksu Kdi; -Kdi -Ksu];    %(matrice totale con ff fb bf bb)
        else
           KOt=ck*(Idelta*depepo+Delta*Kan+Delta1*Kan1+Kost+KTep);
           P=[KOt KOt; -KOt -KOt];
        end

        Mo=diag([Belo; -Belo])+P;

%         'ver exp'
%         keyboard
        Oo=expm_mio(Mo,ifpla,iama,Pust);
%         'ver exp'
%         keyboard

        freqs=freq;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Se la struttura trasversa e' banale   (omogenea)

 else % ai
   if ifp>=1
      display('sono entrata della parte senza struttura trasversa'), pausak
   end


   dos=Li(istr);
   Delta=ani(istr)/rr^2;
   Delta1=ani1(istr)/rr^2;
     if ifp>=1
       disp('[spessore nint larghezza anisotropia]')
       [dos rint ai(istr) Delta]
       pausak
     end

    if dos~=0
       depe=(rint^2-rr^2)/rr^2;
       depez=1-rr^2/rint^2;
       ck=-j*kcav*dos;
       Beli=ck*be;
       Mod=([Beli; -Beli]);
       if iztm==1
          KOt=ck*(Idelta*depe+Delta*Kan+Delta1*Kan1+KTep);
          KOz=ck*Ideltaz*depez+KTepz;
          Ksu=KOt+KOz;
          Kdi=KOt-KOz;
          P=[Ksu Kdi; -Kdi -Ksu];
       else
          KOt=ck*(Idelta*depe+Delta*Kan+Delta1*Kan1+KTep);
          P=[KOt KOt; -KOt -KOt];
       end

       if ifp>=1
         display('matrice K'), pausak
         figure,  surf(abs(KOt)), shading('interp'), view(0,90), colorbar, pausak
       end

            Mo=diag(Mod)+P;
%'expm_mio 2'
%pausak
            iama=[];
            Oo=expm_mio(Mo,ifpla,iama,Pust);
            freqs=freq;
   else
    Oo=1;
   end

  end  %ai


     if imem==0
      Tc=Oo*Tc;
     else
      Tc=prodmat(Oo,Tc,ifpla,Pust);
     end

  if ifp>0
   figure, surf(abs(Oo)), shading('interp'), view(0,90), colorbar, pausak
  end

 else  %ilayfast

  icousa=icousav(istr);
  Tc=Tstof{icousa};
%  'uso dati vecchi '
%  pausak

 end  %ilayfast

else
  Oo=1;

end  %layskip


if layskip(istr)==0 & ilaymem(istr)==1
   if pol==pvet(1)
    icoustor=icoustor+1;
    icousav(istr)=icoustor;
   else
    icoustor=icousav(istr);
   end
%   icoustor
%   pausak
   if length(Tc)==1
    Tstor(:,:,icoustor)=IdeOo;
   else
    Tstor(:,:,icoustor)=Tc;
   end
   if iproga==1
     ' qui per gamma', keyboard

     Ga1=diag(Gas);
     s=size(Tc);
     l1=s(1)/2;
     l2=s(1)/2+1;
     l3=s(1);
     T11=Tc(1:l1,1:l1);
     T12=Tc(1:l1,l2:l3);
     T21=Tc(l2:l3,1:l1);
     T22=Tc(l2:l3,l2:l3);

     GAU=(T21+T22*Ga1)/(T11+T12*Ga1);
     ' qui per gamma dopo', keyboard


   end

%   [icoustor istr]
%   disp(' emme_mix fin: Tstor'), keyboard
   ifplatot=1*isem;
   T=1;
   Tc=1;
end

[Tfas,Ffas]=eltime(Tfas,Ffas,istr);

%if aitr~=0
% disp(' emmean '), keyboard
%end
% [ifpla ifplatot]
% toc
% disp('  '),
%disp(' emme_mix '), keyboard
