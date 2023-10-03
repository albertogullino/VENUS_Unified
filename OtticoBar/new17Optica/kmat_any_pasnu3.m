%' kmat_any', keyboard

if mm~=0
 'caso tutti i modi pasnu=1: devo iniziare da mm=0'
 if ifp==-10
 % keyboard
 end 
end
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



    if iany==1 | ianys==1 | ifirstgrat==1
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

        isub=0;
        inmod=1+(nubesu-nubesi)/pasnu;
        inmodi=1;

  
       
%'qui', keyboard       

%%       if nupd~=0  %(m dispari 1,3,5...)

          Ideldiae=[ Ipeze/ddu; Ideltadiag];
          Ideldiam=[ Ipezm/ddu; Ideltadiag];
%          Ideldiae=[ Ipeze; Ideltadiag];
%          Ideldiam=[ Ipezm; Ideltadiag];
          Ideltea=[ Ipeze; Idelteac];
          Ideltma=[ Ipezm; Ideltmac];

        % even
          Idelteep{1}=(diag(Idelteah,sm1)+diag(Ideltea,-sm1)+diag(Ideldiae))/4;
          Ideltmmp{1}=(-(diag(Ideltmah,sm1)+diag(Ideltma,-sm1))+diag(Ideldiam))/4;
          Ideltemp{1}=-segem*(-diag(Idelteah,sm1)+diag(Ideltea,-sm1)-diag(Ideldiae))/4;
          Ideltmep{1}=-segem*(diag(Ideltmah,sm1)-diag(Ideltma,-sm1)-diag(Ideldiam))/4;

%          Ideltee=Ideltee(iip,iip);
%          Ideltmm=Ideltmm(iip,iip);
%          Ideltem=segem*Ideltem(iip,iip);
%          Ideltme=segem*Ideltme(iip,iip);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Matrice tale che Kt=Kt+Delta(strato)*Ideltaap
%% Cosi' anche le altre per i casi odd e con m pari

%          Ideltaap=[Ideltee -Ideltem; -Ideltme -Ideltmm]/4;
%          Ideltaap1=[Ideltee -Ideltem; -Ideltme Ideltmm]/4;
          Ideltaap1=[Idelteep{1} Ideltemp{1}; Ideltmep{1} Ideltmmp{1}];
%          Ideltaap=[-Ideltee -Ideltem; -Ideltme -Ideltmm]/4;
%          'passo', keyboard
%          Ideltaap=[Ideltee -Ideltem; -Ideltme Ideltmm]/2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % odd
          Idelteem{1}=(diag(Idelteah,sm1)+diag(Ideltea,-sm1)-diag(Ideldiae))/4;
          Ideltmmm{1}=(-(diag(Ideltmah,sm1)+diag(Ideltma,-sm1))-diag(Ideldiam))/4;
          Ideltemm{1}=-segem*(-diag(Idelteah,sm1)+diag(Ideltea,-sm1)+diag(Ideldiae))/4;
          Ideltmem{1}=-segem*(diag(Ideltmah,sm1)-diag(Ideltma,-sm1)+diag(Ideldiam))/4;

%          Ideltee=Ideltee(iip,iip);
%          Ideltmm=Ideltmm(iip,iip);
%          Ideltem=segem*Ideltem(iip,iip);
%          Ideltme=segem*Ideltme(iip,iip);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         PIdeltaam=[Ideltee -Ideltem; -Ideltme Ideltmm]/4;
          Ideltaam1=[Idelteem{1} Ideltemm{1}; Ideltmem{1} Ideltmmm{1}];
%          Ideltaam=[-Ideltee -Ideltem; -Ideltme -Ideltmm]/4;
%          Ideltaam=[Ideltee -Ideltem; -Ideltme Ideltmm]/2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%       else  %nube pari  (m=0,2,4...) !!! Ancora da controllare i conti

       % odd
         fa2=2;
         Idelteauv=[ fa2*Ipeze; Idelteac];       
         Ideltmauv=[ 0*Ipezm; Ideltmac];
         Idelteadv=[ fa2*Ipeze; Idelteac];
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

         Idelteem{2}=(diag(Idelteau,sm1)+diag(Idelteadv,-sm1))/4;
         Ideltmmm{2}=(-(diag(Ideltmau,sm1)+diag(Ideltmadv,-sm1)))/4;
         Ideltemm{2}=segem*(diag(Idelteau,sm1)-diag(Ideltmadv,-sm1))/4;
         Ideltmem{2}=segem*(-diag(Ideltmau,sm1)+diag(Idelteadv,-sm1))/4;
%'qui fat2', keyboard
%         Ideltee=diag(Idelteau0,sm1)+diag(Idelteadv0,-sm1);
%         Ideltmm=-(diag(Ideltmau0,sm1)+diag(Ideltmadv0,-sm1));
%         Ideltem=diag(Idelteau0,sm1)-diag(Ideltmadv0,-sm1);
%         Ideltme=-diag(Ideltmau0,sm1)+diag(Idelteadv0,-sm1);

%          Ideltee=Ideltee(iip,iip);
%          Ideltmm=Ideltmm(iip,iip);
%          Ideltem=segem*Ideltem(iip,iip);
%          Ideltme=segem*Ideltme(iip,iip);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         Ideltaam0=[Idelteem{2} Ideltemm{2}; Ideltmem{2} Ideltmmm{2}];
%         Ideltaam=[Ideltee Ideltem; Ideltme Ideltmm]/2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

       % even

         Idelteauv=[ 0*Ipeze; Idelteac];
         Ideltmauv=[ fa2*Ipezm; Ideltmac];
         Idelteadv=[ 0*Ipeze; Idelteac];
         Ideltmadv=[ fa2*Ipezm; Ideltmac];

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

         Idelteep{2}=(diag(Idelteau,sm1)+diag(Idelteadv,-sm1))/4;
         Ideltmmp{2}=(-(diag(Ideltmau,sm1)+diag(Ideltmadv,-sm1)))/4;
         Ideltemp{2}=segem*(+diag(Idelteau,sm1)-diag(Ideltmadv,-sm1))/4;
         Ideltmep{2}=segem*(-diag(Ideltmau,sm1)+diag(Idelteadv,-sm1))/4;

%         Ideltee=diag(Idelteau0,sm1)+diag(Idelteadv0,-sm1);
%         Ideltmm=-(diag(Ideltmau0,sm1)+diag(Ideltmadv0,-sm1));
%         Ideltem=+diag(Idelteau0,sm1)-diag(Ideltmadv0,-sm1);
%         Ideltme=-diag(Ideltmau0,sm1)+diag(Idelteadv0,-sm1);

%          Ideltee=Ideltee(iip,iip);
%          Ideltmm=Ideltmm(iip,iip);
%          Ideltem=segem*Ideltem(iip,iip);
%          Ideltme=segem*Ideltme(iip,iip);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         PIdeltaap=[Ideltee Ideltem; Ideltme Ideltmm]/4;
         Ideltaap0=[Idelteep{2} Ideltemp{2}; Ideltmep{2} Ideltmmp{2}];
%         Ideltaap=[Ideltee Ideltem; Ideltme Ideltmm]/2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%       end

% ricostruisco le matrici
         iip=1:lk;
         ze=zeros(lk,lk);
         iip4=1:4*lk;
         
  % parte p

         Iee0=Idelteep{2};
         Imm0=Ideltmmp{2};
         Iem0=Ideltemp{2};
         Ime0=Ideltmep{2};
         Iee1=Idelteep{1};
         Imm1=Ideltmmp{1};
         Iem1=Ideltemp{1};
         Ime1=Ideltmep{1};         
         
clear IEE IMM IME IEM
         for imobx=0:1:numodiacc-1
          pux=imobx*lk*2+iip4;
          iis=iip+lk*imobx;
          iis1=iip+lk*(imobx+1);
          
          Iee01=[Iee0(iis,iis) ze Iee0(iis,iis1) ze; ze Iee1(iis,iis) ze Iee1(iis,iis1); ...
          Iee0(iis1,iis) ze Iee0(iis1,iis1) ze; ze Iee1(iis1,iis) ze Iee1(iis1,iis1)];
          
          Imm01=[Imm0(iis,iis) ze Imm0(iis,iis1) ze; ze Imm1(iis,iis) ze Imm1(iis,iis1); ...
          Imm0(iis1,iis) ze Imm0(iis1,iis1) ze; ze Imm1(iis1,iis) ze Imm1(iis1,iis1)];
          
          Iem01=[Iem0(iis,iis) ze Iem0(iis,iis1) ze; ze Iem1(iis,iis) ze Iem1(iis,iis1); ...
          Iem0(iis1,iis) ze Iem0(iis1,iis1) ze; ze Iem1(iis1,iis) ze Iem1(iis1,iis1)];
          
          Ime01=[Ime0(iis,iis) ze Ime0(iis,iis1) ze; ze Ime1(iis,iis) ze Ime1(iis,iis1); ...
          Ime0(iis1,iis) ze Ime0(iis1,iis1) ze; ze Ime1(iis1,iis) ze Ime1(iis1,iis1)];        


          IEE(pux,pux)=Iee01;
          IMM(pux,pux)=Imm01;
          IEM(pux,pux)=Iem01;
          IME(pux,pux)=Ime01;
          
%          pux1=(imobx+2)*lk*2+iip4;
%          iis=iip+lk*imobx;
%          iis1=iip+lk*(imobx+1);
          
%          Iee01=[Iee0(iis,iis) ze Iee0(iis,iis1) ze; ze Iee1(iis,iis) ze Iee1(iis,iis1); ...
%          Iee0(iis1,iis) ze Iee0(iis1,iis1) ze; ze Iee1(iis1,iis) ze Iee1(iis1,iis1)];
%          
%          Imm01=[Imm0(iis,iis) ze Imm0(iis,iis1) ze; ze Imm1(iis,iis) ze Imm1(iis,iis1); ...
%          Imm0(iis1,iis) ze Imm0(iis1,iis1) ze; ze Imm1(iis1,iis) ze Imm1(iis1,iis1)];
%          
%          Iem01=[Iem0(iis,iis) ze Iem0(iis,iis1) ze; ze Iem1(iis,iis) ze Iem1(iis,iis1); ...
%          Iem0(iis1,iis) ze Iem0(iis1,iis1) ze; ze Iem1(iis1,iis) ze Iem1(iis1,iis1)];
%          
%          Ime01=[Ime0(iis,iis) ze Ime0(iis,iis1) ze; ze Ime1(iis,iis) ze Ime1(iis,iis1); ...
%          Ime0(iis1,iis) ze Ime0(iis1,iis1) ze; ze Ime1(iis1,iis) ze Ime1(iis1,iis1)];                  
          
          
          
          
%imobx, pausak          
         end
         pr=1:(numodiacc+1)*lk;
         Ideltaap=[IEE(pr,pr) IEM(pr,pr); IME(pr,pr) IMM(pr,pr)];
         
  % parte m

         Iee0=Idelteem{2};
         Imm0=Ideltmmm{2};
         Iem0=Ideltemm{2};
         Ime0=Ideltmem{2};
         Iee1=Idelteem{1};
         Imm1=Ideltmmm{1};
         Iem1=Ideltemm{1};
         Ime1=Ideltmem{1};         

clear IEE IMM IME IEM

         for imobx=0:1:numodiacc-1
          pux=imobx*lk*2+iip4;
          iis=iip+lk*imobx;
          iis1=iip+lk*(imobx+1);

          Iee01=[Iee0(iis,iis) ze Iee0(iis,iis1) ze; ze Iee1(iis,iis) ze Iee1(iis,iis1); ...
          Iee0(iis1,iis) ze Iee0(iis1,iis1) ze; ze Iee1(iis1,iis) ze Iee1(iis1,iis1)];
          
          Imm01=[Imm0(iis,iis) ze Imm0(iis,iis1) ze; ze Imm1(iis,iis) ze Imm1(iis,iis1); ...
          Imm0(iis1,iis) ze Imm0(iis1,iis1) ze; ze Imm1(iis1,iis) ze Imm1(iis1,iis1)];
          
          Iem01=[Iem0(iis,iis) ze Iem0(iis,iis1) ze; ze Iem1(iis,iis) ze Iem1(iis,iis1); ...
          Iem0(iis1,iis) ze Iem0(iis1,iis1) ze; ze Iem1(iis1,iis) ze Iem1(iis1,iis1)];
          
          Ime01=[Ime0(iis,iis) ze Ime0(iis,iis1) ze; ze Ime1(iis,iis) ze Ime1(iis,iis1); ...
          Ime0(iis1,iis) ze Ime0(iis1,iis1) ze; ze Ime1(iis1,iis) ze Ime1(iis1,iis1)];          
          
          IEE(pux,pux)=Iee01;
          IMM(pux,pux)=Imm01;
          IEM(pux,pux)=Iem01;
          IME(pux,pux)=Ime01;
         end
         Ideltaam=[IEE(pr,pr) IEM(pr,pr); IME(pr,pr) IMM(pr,pr)];
                  Iaam=Ideltaam;
                  Iaap=Ideltaap;
      if ipolar==0
       II=[Ideltaap Ideltaap*0; Ideltaap*0 Ideltaam];
       Ideltaap=II;
       Ideltaam=II;
       Ideltad=[Ideltad Ideltad];
       Ideltazd=[Ideltazd Ideltazd];
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