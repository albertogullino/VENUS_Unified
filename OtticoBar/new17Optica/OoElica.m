nScalini=nScal-1;

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
%              KrV=KOS(:,:,ipri)*depepf;
              Kvi(:,:,ipri)=KOS(:,:,ipri);
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
                 Krz=dpes1*(KOS(:,:,ipri)-KOS(:,:,ipriv))*dpes2*depepfz;
                else
                 Krz=dpes1*KOS(:,:,ipri)*dpes2*depepfz;
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
           KOt=ck*(IDE+Delta*Kan+Delta_gr*Kan_gr+Delta1*Kan1+Kost+KTep);
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
           KOt=ck*(IDE+Delta_gr*Kan_gr+Delta*Kan+Delta1*Kan1+Kost+KTep);
           P=[KOt KOt; -KOt -KOt];
        end
       if ifp==1
        [istr; aitr], 'emme_mix', pausak
       end
 if istocont==1
       ' Emme P', keyboard
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