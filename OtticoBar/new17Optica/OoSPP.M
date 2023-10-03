%if Ps.tuttefreq==0
'dentro OoSPP', keyboard
NsMet=0;
if isfield(Par,'NsMet')==1
 NsMet=Par.NsMet;
 NsNomet=nScalini-NsMet;
end


if shai(1,istr+1)==11
 Puscalini=2:NsNomet;
else
 if NsMet>0
  Puscalini=NsNomet+1:nScalini;
 else 
  Puscalini=2:nScalini;
 end 
end
Kost0=Kost;

if Ps.tuttefreq==1 | ifr==1

          KOS=Kos{sh0};
          KOZ=Kosz{sh0};

          ck=-j*kcav*dos;
          oov{1}=Oo;
%'anota', keyboard


%for iscal=2:nScalini   % in primo scalino è fatto nel main
for iscal=Puscalini   % in primo scalino è fatto nel main


           Kost=Kost0+dpes1*KOS(:,:,iscal)*dpes2*depepf;

         if iztm==1
            KOt=ck*(IDE+ianisotropia*(Delta*Kan+Delta_gr*Kan_gr+Delta1*Kan1)+Kost+KTep);
            Kostz=dpes1*KOZ(:,:,iscal)*dpes2*depepfz;
            
            KOz=ck*(IDEz+Kostz+KTepz);
            P=[KOt+KOz  KOt-KOz; -(KOt-KOz)  -(KOt+KOz)];
        else
           KOt=ck*(IDE+ianisotropia*(Delta_gr*Kan_gr+Delta*Kan+Delta1*Kan1)+Kost+KTep);
           P=[KOt KOt; -KOt -KOt];
        end


% '        passo qui OoSPP', keyboard
      	  P=diag(Belod)+P;

             Ood=expm(P);
             oov{iscal}=Ood;
             Oo=Ood*Oo;
             
end       
 'Calcolata Spirale '
 if pol==1
  OoScala1=Oo;
 else
  OoScala2=Oo;
 end
 OoScala=Oo;
 saK
%save gauss
%save gaussOo OoScala
%if ifp==-10
% 'fine  OoSPP', keyboard
%end 
else  % leggo i dati
 if pol==1
  Oo=OoScala1;
 else
  Oo=OoScala2;
 end

end


 'fine  OoSPP', keyboard
