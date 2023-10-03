
          sd=size(ipar);
          dua=zeros(length(dilu)-1,sd(2),sd(3));
          duze=zeros(size(dilu));
          dua2=zeros(length(dilu2)-1,sd(2),sd(3));
          duze2=zeros(size(dilu2));
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

%'lensas 2 cont',keyboard
          dv=[dv(apu_p); dilu; dv(apu_d); dilu2; dv(bpu_d)];
          ipar=[ipar(apu_p1,:,:); dua; ipar(apu_d,:,:); dua2; ipar(bpu_d,:,:)];
          anyf=[anyf(apu_p); duze; anyf(apu_d); duze2; anyf(bpu_d)];
          fst=[fst(apu_p,:); pu_fst; fst(apu_d,:); pu_fst2; fst(bpu_d,:)];
          iauto=[iauto(apu_p,:); [duze duze]; iauto(apu_d,:); [duze2 duze2]; iauto(bpu_d,:) ];
          ifield=[ifield(apu_p,:); duze; ifield(apu_d,:); duze2; ifield(bpu_d,:)];

          pudum=ones(size(dilu(1:end-1)))*(apu_p1(end)+1);
          pudum2=ones(size(dilu2(1:end-1)))*(apu_p1(end)+1);
          radis=radii.a;

   %       if Ndisc>1 
          if Ndisc>1 | Ndisc<0
           sail=size(ailu);
           sa=1:sail(2);
           said=size(nilu);
           san=1:said(2);
           sail2=size(ailu2);
           if sail2(2)<sail(2)
            ailu2=[ailu2 zeros(sail2(1),sail(2)-sail2(2))];
            nilu2=[nilu2 zeros(sail2(1),sail(2)-sail2(2))];
            sail2=size(ailu2);
           end
           sa2=1:sail(2);
           said2=size(nilu2);
           san2=1:said(2);
           radii.a=[radii.a(apu_p,sa); ailu; radii.a(apu_d,sa); ailu2; radii.a(bpu_d,sa2)];
           radii.b=[radii.b(apu_p,sa); bilu; radii.b(apu_d,sa); bilu2; radii.b(bpu_d,sa2)];
           radii.c=[radii.c(apu_p,sa); ailu*0; radii.c(apu_d,sa); ailu2*0; radii.c(bpu_d,sa2)];
           shavet=[shavet(apu_p,sa); ones(size(ailu))*ShMv(1); shavet(apu_d,sa); ones(size(ailu2))*ShMv(2); shavet(bpu_d,sa)];
           xm=[xm(apu_p,san); -10*ones(size(nilu)); xm(apu_d,san); -10*ones(size(nilu2)); xm(bpu_d,san2)];
           radii.array=[radii.array(apu_p1,:); radii.array(pudum,:); radii.array(apu_d,:); radii.array(pudum2,:); radii.array(bpu_d,:)];
           nv0=[nv0(apu_p,san); nilu; nv0(apu_d,san); nilu2; nv0(bpu_d,san2)];
          else
           radii.a=[radii.a(apu_p,1); ailu; radii.a(apu_d,1); ailu2; radii.a(bpu_d,1)];
           radii.b=[radii.b(apu_p,1); bilu; radii.b(apu_d,1); ailu2*0; radii.b(bpu_d,1)];
           radii.c=[radii.c(apu_p,1); ailu*0; radii.c(apu_d,1); ailu2*0; radii.c(bpu_d,1)];
           shavet=[shavet(apu_p,1); ones(size(ailu))*0; shavet(apu_d,1); ones(size(ailu2))*0; shavet(bpu_d,1)];
%          ' dopo  reass', keyboard
           xm=[xm(apu_p,1); -10*ones(size(nilu)); xm(apu_d,1); -10*ones(size(nilu2)); xm(bpu_d,1)];
           radii.array=[radii.array(apu_p1,1); radii.array(pudum,1); radii.array(apu_d,1); radii.array(pudum2,1); radii.array(bpu_d,1)];
           nva=repmat(nilu,1,2);
           nva(:,2)=0;
           nva2=repmat(nilu2,1,2);
           nva2(:,2)=0;
           nv0=[nv0(apu_p,1:2); nva; nv0(apu_d,1:2); nva2; nv0(bpu_d,1:2)];
          end
shavet_le=shavet(:,1);

%          ' dopo  reass', keyboard
%          ' qui sub', pausak
