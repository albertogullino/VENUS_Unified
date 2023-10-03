%' kmat_any', keyboard
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%/%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Contributo dell'Anisotropia  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('ifirstgrat')
 ifirstgrat=0;
end


 Ideltaap=0;
 Ideltaam=0;
 Ideltaa=0;
 istatt=0;
 flp=1-iLP;
 mr=1;
 Ide=diag(ones(1,2*length(be)));

 if iLP==0
    bb=conj(sqrt(1-flp*KKv.^2*mr^2));
    ZEv=real(mr./bb);
    ZMv=real(mr.*bb);
    ZEv=(mr./bb);
    ZMv=(mr.*bb);
    Ideltad=([ZEv; ZMv])/2;
    Ideltazd=([ZEv*0; ZMv.*KKv.^2./(1-KKv.^2)])/2;
    KKtz=[KKv*0; KKv];
    KKt=[KKv; KKv];
    KKth=[KKv; KKv*0];
    if ipolar==0
     KKtz=[KKtz; KKtz];
     KKth=[KKth; KKth];
     KKt=[KKt; KKt];
     Ideltad=[Ideltad; Ideltad];
     Ideltazd=[Ideltazd; Ideltazd];
    end
    KKtz=KKtz(Pusas);
    KKth=KKth(Pusas);
%    KKtz=KKtz(Pus);
%    KKth=KKth(Pus);


    if iany>=1 | ianys==1 | ifirstgrat==1
       ZEva=mr./sqrt(1-flp*KKva.^2*mr^2);
       ZMva=mr.*sqrt(1-flp*KKva.^2*mr^2);
       Idelteac=ZEva/2;
       Ideltmac=ZMva/2;

       sm1=length(KK);
       Ze=mr./sqrt(1-flp*KK.^2*mr^2);
       Zm=mr.*sqrt(1-flp*KK.^2*mr^2);
       Ipeze=Ze/2;
       Ipezm=Zm/2;

       Ideltea=[ Ipeze; Idelteac];
       Ideltma=[ Ipezm; Ideltmac];
       ddu=1;
       Idelteah=[ Ipeze/ddu; Idelteac];
       Ideltmah=[ Ipezm/ddu; Ideltmac];
       Ideltadiag=zeros(size(Ideltea));

       nupd=nube/2-fix(nube/2);
       if length(mm)==1
        isub=0;
        inmod=1+(nubesu-nubesi)/pasnu;
        inmodi=1;
       else
        isub=1;
        inmod=1+(nubesu-nubesi)/pasnu;
%        inmodi=(fix((nube)/2)-numodiacc)*length(KK)+1;
        inmodi=(fix((nubesi)/2))*length(KK)+1;
       end
%       iipv=[(numodi-2*numodiacc-isub)*length(KK)+1:numodi*length(KK)];
%       keyboard

       iip=[inmodi:inmod*length(KK)];

       if nupd~=0  %(m dispari 1,3,5...)

          Ideldiae=[ Ipeze/ddu; Ideltadiag];
          Ideldiam=[ Ipezm/ddu; Ideltadiag];
%          Ideldiae=[ Ipeze; Ideltadiag];
%          Ideldiam=[ Ipezm; Ideltadiag];
          Ideltea=[ Ipeze; Idelteac];
          Ideltma=[ Ipezm; Ideltmac];

        % even
          Ideltee=diag(Idelteah,sm1)+diag(Ideltea,-sm1)+diag(Ideldiae);
          Ideltmm=-(diag(Ideltmah,sm1)+diag(Ideltma,-sm1))+diag(Ideldiam);
          Ideltem=-diag(Idelteah,sm1)+diag(Ideltea,-sm1)-diag(Ideldiae);
          Ideltme=diag(Ideltmah,sm1)-diag(Ideltma,-sm1)-diag(Ideldiam);

          Ideltee=Ideltee(iip,iip);
          Ideltmm=Ideltmm(iip,iip);
          Ideltem=segem*Ideltem(iip,iip);
          Ideltme=segem*Ideltme(iip,iip);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Matrice tale che Kt=Kt+Delta(strato)*Ideltaap
%% Cosi' anche le altre per i casi odd e con m pari

%          Ideltaap=[Ideltee -Ideltem; -Ideltme -Ideltmm]/4;
          Ideltaap=[Ideltee -Ideltem; -Ideltme Ideltmm]/4;
%          Ideltaap=[-Ideltee -Ideltem; -Ideltme -Ideltmm]/4;
%          'passo', keyboard
%          Ideltaap=[Ideltee -Ideltem; -Ideltme Ideltmm]/2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % odd
          Ideltee=diag(Idelteah,sm1)+diag(Ideltea,-sm1)-diag(Ideldiae);
          Ideltmm=-(diag(Ideltmah,sm1)+diag(Ideltma,-sm1))-diag(Ideldiam);
          Ideltem=-diag(Idelteah,sm1)+diag(Ideltea,-sm1)+diag(Ideldiae);
          Ideltme=diag(Ideltmah,sm1)-diag(Ideltma,-sm1)+diag(Ideldiam);

          Ideltee=Ideltee(iip,iip);
          Ideltmm=Ideltmm(iip,iip);
          Ideltem=segem*Ideltem(iip,iip);
          Ideltme=segem*Ideltme(iip,iip);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         PIdeltaam=[Ideltee -Ideltem; -Ideltme Ideltmm]/4;
          Ideltaam=[Ideltee -Ideltem; -Ideltme Ideltmm]/4;
%          Ideltaam=[-Ideltee -Ideltem; -Ideltme -Ideltmm]/4;
%          Ideltaam=[Ideltee -Ideltem; -Ideltme Ideltmm]/2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

       else  %nube pari  (m=0,2,4...) !!! Ancora da controllare i conti

       % odd
         Idelteauv=[ 2*Ipeze; Idelteac];       
         Ideltmauv=[ 0*Ipezm; Ideltmac];
         Idelteadv=[ 2*Ipeze; Idelteac];
         Ideltmadv=[ 0*Ipezm; Ideltmac];

         Idelteau=[ 1*Ipeze; Idelteac];
         Ideltmau=[ 0*Ipezm; Ideltmac];
         Ideltead=[ 1*Ipeze; Idelteac];
         Ideltmad=[ 0*Ipezm; Ideltmac];
         Idp=  [ Ipeze; Idelteac;   0*Ipezm; Ideltmac];
         Idm=  [ 0*Ipeze; Idelteac;   Ipezm; Ideltmac];

         Idelteauv0=[ 0*Ipeze; Idelteac];
         Ideltmauv0=[ 0*Ipezm; Ideltmac];
         Idelteadv0=[ 0*Ipeze; Idelteac];
         Ideltmadv0=[ 0*Ipezm; Ideltmac];

         Idelteau0=[ 0*Ipeze; Idelteac];
         Ideltmau0=[ 0*Ipezm; Ideltmac];
         Ideltead0=[ 0*Ipeze; Idelteac];
         Ideltmad0=[ 0*Ipezm; Ideltmac];

         Ideltee=diag(Idelteau,sm1)+diag(Idelteadv,-sm1);
         Ideltmm=-(diag(Ideltmau,sm1)+diag(Ideltmadv,-sm1));
         Ideltem=diag(Idelteau,sm1)-diag(Ideltmadv,-sm1);
         Ideltme=-diag(Ideltmau,sm1)+diag(Idelteadv,-sm1);

%         Ideltee=diag(Idelteau0,sm1)+diag(Idelteadv0,-sm1);
%         Ideltmm=-(diag(Ideltmau0,sm1)+diag(Ideltmadv0,-sm1));
%         Ideltem=diag(Idelteau0,sm1)-diag(Ideltmadv0,-sm1);
%         Ideltme=-diag(Ideltmau0,sm1)+diag(Idelteadv0,-sm1);

          Ideltee=Ideltee(iip,iip);
          Ideltmm=Ideltmm(iip,iip);
          Ideltem=segem*Ideltem(iip,iip);
          Ideltme=segem*Ideltme(iip,iip);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         Ideltaam=[Ideltee Ideltem; Ideltme Ideltmm]/4;
%         Ideltaam=[Ideltee Ideltem; Ideltme Ideltmm]/2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

       % even
         Idelteauv=[ 0*Ipeze; Idelteac];
         Ideltmauv=[ 2*Ipezm; Ideltmac];
         Idelteadv=[ 0*Ipeze; Idelteac];
         Ideltmadv=[ 2*Ipezm; Ideltmac];

         Idelteau=[ 0*Ipeze; Idelteac];
         Ideltmau=[ 1*Ipezm; Ideltmac];
         Ideltead=[ 0*Ipeze; Idelteac];
         Ideltmad=[ 1*Ipezm; Ideltmac];

         Idelteauv0=[ 0*Ipeze; Idelteac];
         Ideltmauv0=[ 0*Ipezm; Ideltmac];
         Idelteadv0=[ 0*Ipeze; Idelteac];
         Ideltmadv0=[ 0*Ipezm; Ideltmac];

         Idelteau0=[ 0*Ipeze; Idelteac];
         Ideltmau0=[ 0*Ipezm; Ideltmac];
         Ideltead0=[ 0*Ipeze; Idelteac];
         Ideltmad0=[ 0*Ipezm; Ideltmac];

         Ideltee=diag(Idelteau,sm1)+diag(Idelteadv,-sm1);
         Ideltmm=-(diag(Ideltmau,sm1)+diag(Ideltmadv,-sm1));
         Ideltem=+diag(Idelteau,sm1)-diag(Ideltmadv,-sm1);
         Ideltme=-diag(Ideltmau,sm1)+diag(Idelteadv,-sm1);

%         Ideltee=diag(Idelteau0,sm1)+diag(Idelteadv0,-sm1);
%         Ideltmm=-(diag(Ideltmau0,sm1)+diag(Ideltmadv0,-sm1));
%         Ideltem=+diag(Idelteau0,sm1)-diag(Ideltmadv0,-sm1);
%         Ideltme=-diag(Ideltmau0,sm1)+diag(Idelteadv0,-sm1);

          Ideltee=Ideltee(iip,iip);
          Ideltmm=Ideltmm(iip,iip);
          Ideltem=segem*Ideltem(iip,iip);
          Ideltme=segem*Ideltme(iip,iip);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         PIdeltaap=[Ideltee Ideltem; Ideltme Ideltmm]/4;
         Ideltaap=[Ideltee Ideltem; Ideltme Ideltmm]/4;
%         Ideltaap=[Ideltee Ideltem; Ideltme Ideltmm]/2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       end
       
       if iany_CIRC==1
       
          Ideldiae=diag(Ideltea);
          Ideldiam=diag(Ideltma);
          IdeltaaC=[Ideldiae 0*Ideldiae; 0*Ideldiam Ideldiam];          
       
       end
       
       
%' any ' ,keyboard
%
%         Ideltaaps=Ideltaap;
%         Ideltaams=Ideltaam;

         Ideltaap=Ideltaap(Pusas,Pusas);
         Ideltaam=Ideltaam(Pusas,Pusas);
  end %iany==1
 Ideltad2=Ideltad;
 Ideltazd2=Ideltazd;
 else  %iLP~=0
    ZEv=ones(size(KKv));
    Ideltad2=ZEv2/2;
    Ideltad=ZEv/2;
    Ideltazd=([ZEv.*KKv.^2./(1-KKv.^2)])/2;
    Ideltazd2=([ZEv2.*KKvd.^2./(1-KKvd.^2)])/2;
    KKt=KKv;
 end
%  if ipolar==0
%   Ideltad=[Ideltad; Ideltad];
%   Ideltazd=[Ideltazd; Ideltazd];
%  end
% Matd=Ideltad./pes;
% Matdz=Ideltazd./pes;
 Matd=Ideltad2;
 Matdz=Ideltazd2;
 Mixmat=zeros(size(Matd));
 Mixmatz=zeros(size(Matd));
% ' mat;',keyboard
 if irid_bas==1
  if exist('Pusasf')
   Mixmat(Pusasf)=Matd(Pusasf);
   Mixmatz(Pusasf)=Matdz(Pusasf);
   Madd0=diag(Mixmat(Pusas));
   Maddz0=diag(Mixmatz(Pusas));
  else
   Madd0=0;
   Maddz0=0;
  end
  Ideltad=Ideltad2(Pusas);
  Ideltazd=Ideltazd2(Pusas);
 end


 Idelta=diag(Ideltad);
 Ideltaz=diag(Ideltazd);
 
Idelta_p=Idelta; 
Idelta_m=Idelta;
Ideltaz_p=Ideltaz; 
Ideltaz_m=Ideltaz;
Idep=Ideltazd;
if ~exist('nupd')
 nupd=1;
end
if nupd==10
 Idelta_m=diag(Idp);
 Idelta_p=diag(Idm);
 fi=find(Idp==0);
 Idelpp=Idep;
 Idelpp(fi)=0;
 Ideltaz_p=diag(Idelpp); 
 fi=find(Idm==0);
 Idelpm=Idep;
 Idelpm(fi)=0;  
 Ideltaz_m=diag(Idelpm); 
end

% if iLP==1
%  Ideltaz=0;
% end
% keyboard

%kv=KKvd(Pusas);

%if exist('pedk')
% fapek=(1-1./(1+(kv/kl).^exp_los))*perdk*(-j);
%else
% fapek=zeros(size(kv));
%end
fapek=0;

if ifp>=-3, disp(' Kmat_any '), if ifp>1, keyboard, end, end

 if ired_ret>=2
  Par{1}.dpes1=dpes1;
  Par{1}.dpes2=dpes2; 
  Par{1}.Madd0=Madd0; 
  Par{1}.Maddz0=Maddz0; 
  Par{1}.Idelta=Idelta;  
  Par{1}.Ideltaz=Ideltaz; 
  Par{1}.be=be(Pusas);  
   dPe=diag(Par{1}.dpes1); 
   diatut=dPe;
   Pured_ret0=find(abs(diatut-diatut(1))<1e-6);

   diatut=[dPe; dPe];
   Pured_ret=find(abs(diatut-diatut(1))<1e-6);
   Par{2}.Pured_ret0=Pured_ret0; 
   Par{2}.Pured_ret=Pured_ret; 
   
  dup=Par{1}.Pusas;
  du=dup(Pured_ret0);
  dMadd0=diag(Mixmat(du));
  dMaddz0=diag(Mixmatz(du)); 
  Par{2}.Madd0=dMadd0; 
  Par{2}.Maddz0=dMaddz0;  
  dIdelta=Idelta(Pured_ret0,Pured_ret0); 
  Par{2}.Idelta=dIdelta;  
  dIdeltaz=Ideltaz(Pured_ret0,Pured_ret0); 
  Par{2}.Ideltaz=dIdeltaz; 
  Par{2}.be=be(du); 
 dpes=diag(pes(du)); 
 if iant==0
  dpes1=1;
  dpes2=dpes;
 else
  dpes2=1;
  dpes1=dpes;
 end
 Par{2}.dpes1=dpes1;
 Par{2}.dpes2=dpes2;
 KKtz=KKtz(Pured_ret0);
 KKth=KKth(Pured_ret0);
 if length(KKv)>length(Pured_ret0)
  KKv=KKv(Pured_ret0);
 end 
 KK =Par{2}.KK;
 lKAn =Par{2}.lKAn;
 KKt=[KK; KK];
 numodi=Par{2}.numodi;

 end

%'fine Kmat_any', keyboard



if ifp~=-4
%'QUI any ', keyboard
end