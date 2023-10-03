
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

  %        ' sail prima', keyboard

          dv=[dv(pu_p); dilu; dv(pu_d)];
          ipar=[ipar(pu_p1,:,:); dua; ipar(pu_d,:,:)];
          anyf=[anyf(pu_p); duze; anyf(pu_d)];
          fst=[fst(pu_p,:); pu_fst; fst(pu_d,:)];
%          iauto=[iauto(pu_p,:); [zeros(size(dilu)) zeros(size(dilu))]; iauto(pu_d,:) ];

          iau_add=zeros(length(dilu)-1,size(iauto,2));
          iau_add=[iau_add; iauto(pu_d(1)-1,:)];

          iauto=[iauto(pu_p,:); iau_add; iauto(pu_d,:) ];
          ifield=[ifield(pu_p,:); zeros(size(dilu)); ifield(pu_d,:)];

          pudum=ones(size(dilu(1:end-1)))*(pu_p1(end)+1);

%          ' sail lente prima', keyboard
          fifor=length(find(shavet<7 & shavet>0));
          if fifor>0
           ncMin=1;
          end 
          if Ndisc>1 | Ndisc<0
           sail=size(ailu);

%           let=length(pu_p)+length(ailu)+length(pu_d);
%           larg=max([size(ailu,2) size(radii.a,2)]);
%           rapa=zeros(let,larg);
           sai=sail(2);
           sar=size(radii.a,2);
           if sai>sar
            lari=size(radii.a,1);
            sa=1:sail(2);
            said=size(nilu);
            san=1:said(2);
            rapa=[radii.a zeros(lari,sai-sar)];            
            rapb=[radii.b zeros(lari,sai-sar)];            
            rapc=[radii.c zeros(lari,sai-sar)];            
            shapdu=[shavet zeros(lari,sai-sar)];            
            xmp=[xm zeros(lari,sai-sar)];            
            nvd=[nv0 zeros(lari,sai-sar)];
            ailuc=ailu;
            biluc=bilu;
            niluc=nilu;
           else
            sa=1:sar;
            san=1:sar+1;           
            lari=size(ailu,1);
            ailuc=[ailu zeros(lari,sar-sai)];            
            biluc=[bilu zeros(lari,sar-sai)];            
            niluc=[nilu zeros(lari,sar-sai)];            
            rapa=radii.a;            
            rapb=radii.b;            
            rapc=radii.c;            
            shapdu=shavet;            
            xmp=xm;            
            nvd=nv0;
           end
           radii.a=[rapa(pu_p,sa); ailuc; rapa(pu_d,sa)];
           radii.b=[rapb(pu_p,sa); biluc; rapb(pu_d,sa)];
           radii.c=[rapc(pu_p,sa); ailuc*0; rapc(pu_d,sa)];
           shavet=[shapdu(pu_p,sa); ones(size(ailuc))*ShMv(1); shapdu(pu_d,sa)];
           om=-10*ones(size(niluc));
           xm=[xmp(pu_p,san); om(:,san); xmp(pu_d,san)];
           radii.array=[radii.array(pu_p1,:); radii.array(pudum,:); radii.array(pu_d,:)];
           nv0=[nvd(pu_p,san); niluc; nvd(pu_d,san)];

%           radii.a=[radii.a(pu_p,sa); ailu; radii.a(pu_d,sa)];
%           radii.b=[radii.b(pu_p,sa); bilu; radii.b(pu_d,sa)];
%           radii.c=[radii.c(pu_p,sa); ailu*0; radii.c(pu_d,sa)];
%           shavet=[shavet(pu_p,sa); ones(size(ailu))*ShMv(1); shavet(pu_d,sa)];
%           radii.array=[radii.array(pu_p1,:); radii.array(pudum,:); radii.array(pu_d,:)];
%           nv0=[nv0(pu_p,san); nilu; nv0(pu_d,san)];

          else
           radii.a=[radii.a(pu_p,1); ailu; radii.a(pu_d,1)];
           radii.b=[radii.b(pu_p,1); bilu; radii.b(pu_d,1)];
           radii.c=[radii.c(pu_p,1); ailu*0; radii.c(pu_d,1)];
           shavet=[shavet(pu_p,1); ones(size(ailu))*0; shavet(pu_d,1)];
%          ' dopo  reass', keyboard
           xm=[xm(pu_p,1); -10*ones(size(nilu)); xm(pu_d,1)];
           radii.array=[radii.array(pu_p1,1); radii.array(pudum,1); radii.array(pu_d,1)];
           nva=repmat(nilu,1,2);
           nva(:,2)=0;
           nv0=[nv0(pu_p,1:2); nva; nv0(pu_d,1:2)];
          end
shavet_le=shavet(:,1);
%' dopo  lensas1 ', keyboard
%          ' qui sub', pausak
