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
 KKs=KKv;
 fis=find(abs(KKv-1)<.02);
 KKs(fis)=.98;

 if iLP==0
    bb=conj(sqrt(1-flp*KKs.^2*mr^2));
    
    ZEv=real(mr./bb);
    ZMv=real(mr.*bb);
    ZEv=(mr./bb);
    ZMv=mr.*bb;
    Ideltad=([ZEv; ZMv])/2;
    Ideltazd=([ZEv*0; ZMv.*KKv.^2./(1-KKs.^2)])/2;
    KKtz=[KKs*0; KKv];
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


    if iany==1 | ianys==1 | ifirstgrat==1
      KKvas=KKva;
      fis=find(abs(KKva-1)<.02);
      KKvas(fis)=.98;

       ZEva=mr./sqrt(1-flp*KKvas.^2*mr^2);
       ZMva=mr.*sqrt(1-flp*KKva.^2*mr^2);
       Idelteac=ZEva/2;
       Ideltmac=ZMva/2;

       sm1=length(KK);
       Ze=mr./sqrt(1-flp*KKs.^2*mr^2);
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
       if mm<=1
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

       iip=[inmodi:round((nubesu+1)/2)*length(KK)];

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

          Ideltaap=[Ideltee -Ideltem; -Ideltme Ideltmm]/4;
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


%if ifp~=-4
%'QUI any ', keyboard
%end