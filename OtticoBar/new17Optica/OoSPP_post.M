%if Ps.tuttefreq==0
%'dentro OoSPP', keyboard

          KOS=KSpp;
          KOZ=KSppz;


nScalini=7
puscal=[1:7]
%modos=[0.2 0.2 0.4 .8 1.6 2.8 1]*falu; 
%modos=fliplr(modos)
modos=ones(1,7)*falu;
sum(modos)

%         IDE=Idelta*depepo;
%         IDEz=Ideltaz*depepoz;
%                 depepf=(riv(1)^2-riv(2)^2)/rr^2
%                 depepfz=-rr^2/riv(1)^2+rr^2/rext^2   
rext=1;
rs0=1.9;
                 depepf0=(rs0^2-rext^2)/rr^2
                 depepfz0=-rr^2/rs0^2+rr^2/rext^2   
'somma'
%pausak
Oo=1;
for iscal=puscal   % in primo scalino è fatto nel main

           rs=rs0*modos(iscal);
           depepf=(rs^2-rext^2)/rr^2
           depepfz=-rr^2/rs^2+rr^2/rext^2   
           Kost=dpes1*KOS(:,:,iscal)*dpes2*depepf;

         if iztm==1
            KOt=ck*(IDE+Kost);
            Kostz=dpes1*KOZ(:,:,iscal)*dpes2*depepfz;
            
            KOz=ck*(IDEz+Kostz);
            P=[KOt+KOz  KOt-KOz; -(KOt-KOz)  -(KOt+KOz)];
        else
           KOt=ck*(IDE+Kost);
           P=[KOt KOt; -KOt -KOt];
        end


% '        passo qui OoSPP', keyboard
      	  P=diag(Belod)+P;

%             Ood=expm(P*modos(iscal));
             Ood=expm(P);
             oov{iscal}=Ood;
             Oo=Ood*Oo;
             
end       
OoScala=Oo;
% 'Calcolata Spirale ', keyboard
