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
Kost0=Kost;
Kostz0=Kostz;
Oo=1;
if SPP_GA==1
  Puscalini=1:nScalini; 
end
%'dentro OoSPP GA', keyboard


if Ps.tuttefreq==1 | ifr==1

          KOS=Kos{11};
          KOZ=Kosz{11};

          ck=-j*kcav*dos;
'anota', keyboard


%for iscal=2:nScalini   % in primo scalino è fatto nel main

if ScalaARC==1
 Puscalini=[1:ScaInc];
 nARC=Par.narc;
 depepfVa=(nARC^2-rext^2)/rr^2
 depepfzVa=-rr^2/nARC^2+rr^2/rext^2  
end

nrSca=linspace(rext,rint,nScalini);

nrSca=linspace(1.3,1.8,nScalini);
'MODIFICATO INDICE VORTICE!!!!!!!!'

   if iscell(KOS)
    KOS=KOS{1};
    KOZ=KOZ{1};
   end

depepfV=(nrSca.^2-rext^2)/rr^2;
depepfzV=-rr^2./nrSca.^2+rr^2/rext^2;    
Kost=Kost0;
Kostz=Kostz0;
for iscal=Puscalini   % in primo scalino è fatto nel main
         Kostd=dpes1*KOS(:,:,iscal)*dpes2*depepfV(iscal);
         Kost=Kost+Kostd;
         
         Kostzd=dpes1*KOZ(:,:,iscal)*dpes2*depepfzV(iscal);       
         Kostz=Kostz+Kostzd;       
end       
            KOt=ck*(IDE+ianisotropia*(Delta*Kan+Delta_gr*Kan_gr+Delta1*Kan1)+Kost+KTep);
            KOz=ck*(IDEz+Kostz+KTepz);
            P=[KOt+KOz  KOt-KOz; -(KOt-KOz)  -(KOt+KOz)];


% '        passo qui OoSPP', keyboard
       	    P=diag(Belod)+P;

             Oo=expm(P);
% 'Calcolata Spirale GA ', keyboard
 
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
