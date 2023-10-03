if exist('PApa')==1
   ipudpa=find(shavet_le==9);
%       radii.arrayd{1}=Dp;
%       radii.arrayd{2}=Cp;
%       radii.arrayd{3}=Ndispa; 
%       radii.arrayd{4}=Refp;
%       radii.arrayd{5}=Th;
%       radii.arrayd{6}=Refl;
%       radii.arrayd{7}=0;
   Rp=PApa{1};
   cs=PApa{2};
   Np=PApa{3};
   npar=PApa{4};
   dlay=PApa{5};
   nlay=PApa{6};
   [dilu,ailu,nilu,filu]=particle(Rp,Np,cs,dlay,npar,nlay,ifp)
   bilu=ailu;
           itrle=length(ipudpa);
           pu_p=1:ipudpa(1)-1;
           pu_p1=1:ipudpa(1);
           pu_d=ipudpa(itrle)+1:length(dv);
           pu_fst=filu;

          sd=size(ipar);
          dua=zeros(length(dilu)-1,sd(2),sd(3));
          duze=zeros(size(dilu));
          dv_sa=dv;
          ipar_sa=ipar;
          anyf_sa=anyf;
          fst_sa=fst;
          iauto_sa=iauto;
          ifield_sa=ifield;
          radii_sa=radii;
          sha_sa=shavet;
          n_sa=nv0;
          xm_sa=xm;

%          ' sail prima', keyboard

          dv=[dv(pu_p); dilu; dv(pu_d)];
          ipar=[ipar(pu_p1,:,:); dua; ipar(pu_d,:,:)];
          anyf=[anyf(pu_p); duze; anyf(pu_d)];
          fst=[fst(pu_p,:); pu_fst; fst(pu_d,:)];
          iauto=[iauto(pu_p,:); [zeros(size(dilu)) zeros(size(dilu))]; iauto(pu_d,:) ];
          ifield=[ifield(pu_p,:); zeros(size(dilu)); ifield(pu_d,:)];

          pudum=ones(size(dilu(1:end-1)))*(pu_p1(end)+1);


           radii.a=[radii.a(pu_p,1); ailu; radii.a(pu_d,1)];
           radii.b=[radii.b(pu_p,1); bilu; radii.b(pu_d,1)];
           radii.c=[radii.c(pu_p,1); ailu*0; radii.c(pu_d,1)];
           shavet=[shavet(pu_p,1); ones(size(ailu))*0; shavet(pu_d,1)];
%          ' dopo  reass', keyboard
           xm=[xm(pu_p,1); -10*ones(size(ailu)); xm(pu_d,1)];
           radii.array=[radii.array(pu_p1,1); radii.array(pudum,1); radii.array(pu_d,1)];
           nva=nilu;
           nva(:,2)=0;
           nv0=[nv0(pu_p,1:2); nva; nv0(pu_d,1:2)];

          ' dopo  part_ass ', keyboard
%          ' qui sub', pausak

   
   
end   