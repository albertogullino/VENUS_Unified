%if Ps.tuttefreq==0
%'dentro OoSPP', keyboard
NsMet=0;
if isfield(Par,'NsMet')==1
 NsMet=Par.NsMet;
 NsNomet=nScalini-NsMet;
end


if shai(1,istr+1)==11
 Puscalini=1:NsNomet;
else
 if NsMet>0
  Puscalini=NsNomet+1:nScalini;
 else 
  Puscalini=1:nScalini;
 end 
end
Oo=1;
%'dentro OoSPP', keyboard


if Ps.tuttefreq==1 | ifr==1

          KOS=Kos{11};
          KOZ=Kosz{11};

          ck=-j*kcav*dos;
%'anota', keyboard


%for iscal=2:nScalini   % in primo scalino è fatto nel main

if ScalaARC>=1
 Puscalini=[1:ScaInc];
 nARC=Par.narc;
 depepfVa=(nARC^2-rext^2)/rr^2
 depepfzVa=-rr^2/nARC^2+rr^2/rext^2 
end

if ScalaARC==2
 Puscalini=[1:ScaInc+1];
 t_a=lambda/4/nARC;
 t_a=Par.ARC;
% 'VER ARC', keyboard
%'VERIFICA !!!!!!!!!!!!!!!!!'
% t_a=dos
 t_s=dos;
 if t_a>2*t_s
  'CASO NON SUPPORTATO !!!!!!'
  keyboard
 end
 if t_a>t_s
  delta=t_a-t_s;
  dminor=0;
 else
  delta=t_a;
  dminor=1;
 end
 Delta=t_s-delta; 
 ck_d=-j*kcav*delta; 
 ck_D=-j*kcav*Delta; 
 ck_s=-j*kcav*t_s; 
else 
 dminor=0;
end


%'vedi SPP ARC', keyboard
if ScalaARC==2

 for iscal=Puscalini   % in primo scalino è fatto nel main

          puscaD=1:iscal-1;
          puscad2=iscal:iscal+1; 
          puscad1=iscal; 
          
        if dminor==0  %delta maggiore di ZERO
          
           if iscal==Puscalini(end)
            ck_d=ck_s;
            delta=t_s;
            ck_D=0;
            puscad2=iscal;             
           end
           
           if length(puscaD)>0
            SuKOS=sum(KOS(:,:,puscaD),3);
            SuKOSz=sum(KOZ(:,:,puscaD),3);
           else
            SuKOS=0;
            SuKOSz=0;
           end

% prima parte con Delta   
%'qui debug', keyboard
          ck=ck_D;
          rak=Delta/dos;
          if abs(ck)>0 
           puscad=puscad1;
           SuKOS_d=sum(KOS(:,:,puscad),3);
           SuKOSz_d=sum(KOZ(:,:,puscad),3);

           Kost=dpes1*(SuKOS*depepfV+SuKOS_d*depepfVa)*dpes2;
           Kostz=dpes1*(SuKOSz*depepfzV+SuKOSz_d*depepfzVa)*dpes2;

           KOt=ck*(IDE+ianisotropia*(Delta*Kan+Delta_gr*Kan_gr+Delta1*Kan1)+Kost+KTep);
           KOz=ck*(IDEz+Kostz+KTepz);
           P=[KOt+KOz  KOt-KOz; -(KOt-KOz)  -(KOt+KOz)];
       	   P=rak*diag(Belod)+P;
           Ood=expm(P);
             
           Oo=Ood*Oo;         
          end 
          iscal

%'fine Delta', keyboard
% seconda parte con delta           

           ck=ck_d;
           rak=delta/dos;           
           puscad=puscad2;
           SuKOS_d=sum(KOS(:,:,puscad),3);
           SuKOSz_d=sum(KOZ(:,:,puscad),3);

           Kost=dpes1*(SuKOS*depepfV+SuKOS_d*depepfVa)*dpes2;
           Kostz=dpes1*(SuKOSz*depepfzV+SuKOSz_d*depepfzVa)*dpes2;

           KOt=ck*(IDE+ianisotropia*(Delta*Kan+Delta_gr*Kan_gr+Delta1*Kan1)+Kost+KTep);
           KOz=ck*(IDEz+Kostz+KTepz);
           P=[KOt+KOz  KOt-KOz; -(KOt-KOz)  -(KOt+KOz)];
       	   P=rak*diag(Belod)+P;
           Ood=expm(P);
             
           Oo=Ood*Oo;  
%'fine delta', keyboard           




           
       else %dminor=1
        

           if length(puscaD)>0
            SuKOS=sum(KOS(:,:,puscaD),3);
            SuKOSz=sum(KOZ(:,:,puscaD),3);
           else
            SuKOS=0;
            SuKOSz=0;
           end
           
           iscal

% prima parte con Delta           
          ck=ck_D;
          if iscal==1
           ck=0;
          end      
          rak=Delta/dos;
          if length(SuKOS)>1 & abs(ck)>0
           Kost=dpes1*(SuKOS*depepfV)*dpes2;
           Kostz=dpes1*(SuKOSz*depepfzV)*dpes2;

           KOt=ck*(IDE+ianisotropia*(Delta*Kan+Delta_gr*Kan_gr+Delta1*Kan1)+Kost+KTep);
           KOz=ck*(IDEz+Kostz+KTepz);
           P=[KOt+KOz  KOt-KOz; -(KOt-KOz)  -(KOt+KOz)];
       	   P=rak*diag(Belod)+P;
           Ood=expm(P);
             
           Oo=Ood*Oo;         
%'dopo ck_D', keyboard  
          end


% seconda parte con delta           

           ck=ck_d;
           rak=delta/dos;
           puscad=iscal;
           SuKOS_d=KOS(:,:,puscad);
           SuKOSz_d=KOZ(:,:,puscad);

           Kost=dpes1*(SuKOS*depepfV+SuKOS_d*depepfVa)*dpes2;
           Kostz=dpes1*(SuKOSz*depepfzV+SuKOSz_d*depepfzVa)*dpes2;

           KOt=ck*(IDE+ianisotropia*(Delta*Kan+Delta_gr*Kan_gr+Delta1*Kan1)+Kost+KTep);
           KOz=ck*(IDEz+Kostz+KTepz);
           P=[KOt+KOz  KOt-KOz; -(KOt-KOz)  -(KOt+KOz)];
%           P1=[KOt+KOz  KOt-KOz; -(KOt-KOz)  -(KOt+KOz)];
       	   P=rak*diag(Belod)+P;
           Ood=expm(P);
             
           Oo=Ood*Oo;

%'dopo ck_d', keyboard

                     
       end  %dminor
  
             
 end       %iscal

else   %ScalaARC=1

           iscal=1;
           Kost=dpes1*KOS(:,:,iscal)*dpes2*depepfVa;
           Kostz=dpes1*KOZ(:,:,iscal)*dpes2*depepfzVa;       
	   KOt=ck*(IDE+ianisotropia*(Delta*Kan+Delta_gr*Kan_gr+Delta1*Kan1)+Kost+KTep);
	   KOz=ck*(IDEz+Kostz+KTepz);
	   P=[KOt+KOz  KOt-KOz; -(KOt-KOz)  -(KOt+KOz)];
	   P=diag(Belod)+P;
           Ood=expm(P);
           Oo=Ood;
           Os=Ood;
           save sa Os
           'Vecchio Sopra', keyboard

 for iscal=Puscalini   % in primo scalino è fatto nel main

           Kost=dpes1*KOS(:,:,iscal)*dpes2*depepfV;
           Kostz=dpes1*KOZ(:,:,iscal)*dpes2*depepfzV;       

       
           if ScalaARC==1           
            Kost=Kost+dpes1*KOS(:,:,iscal+ScaInc)*dpes2*depepfVa;
            Kostz=Kostz+dpes1*KOZ(:,:,iscal+ScaInc)*dpes2*depepfzVa;              
           end

            KOt=ck*(IDE+ianisotropia*(Delta*Kan+Delta_gr*Kan_gr+Delta1*Kan1)+Kost+KTep);
            KOz=ck*(IDEz+Kostz+KTepz);
            P=[KOt+KOz  KOt-KOz; -(KOt-KOz)  -(KOt+KOz)];


% '        passo qui OoSPP', keyboard
       	    P=diag(Belod)+P;

             Ood=expm(P);
           Os=Ood;
           save sa Os             
             iscal
             'dopo ck_d VECCHIO', keyboard
             
%             iama=[];
%             ifpla=0;
%             iautm=1;
%             Ood1=expm_mioVor(P,ifpla,iama,Pust,PUrid,iautm);   
             
             oov{iscal}=Ood;
             Oo=Ood*Oo;
             
 end

end  %ScalaARC


 'Calcolata Spirale '
 if pol==1
  OoScala1=Oo;
 else
  OoScala2=Oo;
 end
 OoScala=Oo;
 %saK
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


% 'fine  OoSPP', keyboard
