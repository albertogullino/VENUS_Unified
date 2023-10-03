

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
   dlay=PApa{5}/1000;
   nlay=PApa{6};
   shapar=PApa{7};
%' part', keyboard   
   if length(shapar)==0
    shapar=7;
   end 
   [dilp,ailp,nilp,filp]=particle(Rp,Np,cs,dlay,npar,nlay,ifp);
%             ' dopo  part_ass ', keyboard
%             ' dopo  part_ass ', keyboard

   bilp=ailp;
           itrle=length(ipudpa);
           pu_p=1:ipudpa(1)-1;
           pu_p1=1:ipudpa(1);
           pu_d=ipudpa(itrle)+1:length(dv);
           pu_fst=filp;

          sd=size(ipar);
          dua=zeros(length(dilp)-1,sd(2),sd(3));
          duze=zeros(size(dilp));
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

          dv=[dv(pu_p); dilp; dv(pu_d)];
          ipar=[ipar(pu_p1,:,:); dua; ipar(pu_d,:,:)];
          anyf=[anyf(pu_p); duze; anyf(pu_d)];
          fst=[fst(pu_p,:); pu_fst; fst(pu_d,:)];
          iauto=[iauto(pu_p,:); [zeros(size(dilp)) zeros(size(dilp))]; iauto(pu_d,:) ];
          ifield=[ifield(pu_p,:); zeros(size(dilp)); ifield(pu_d,:)];

          pudum=ones(size(dilp(1:end-1)))*(pu_p1(end)+1);

          if Ndisc>1 | Ndisc<0
           sail=size(ailp);
           sa=1:sail(2);
           said=size(nilp);
           san=1:said(2);
           sai=sail(2);
           sar=size(radii.a,2);
           lari=size(radii.a,1);
           if sai>sar
            rapa=[radii.a zeros(lari,sai-sar)];            
            rapb=[radii.b zeros(lari,sai-sar)];            
            rapc=[radii.c zeros(lari,sai-sar)];            
            shapdu=[shavet zeros(lari,sai-sar)];            
            xmp=[xm zeros(lari,sai-sar)];            
            nvd=[nv0 zeros(lari,sai-sar)];
           else
            sa=1:sar;
%            snv=size(nv0);
%            san=1:snv(2);
            san=1:sar+1;
            ailsa=ailp;
            ailp=[ailsa zeros(sail(1),-sai+sar)];            
            nilsa=nilp;
            nilp=[nilsa zeros(sail(1),-sai+sar)];    
            bilp=ailp;
            rapa=radii.a;            
            rapb=radii.b;            
            rapc=radii.c;            
            shapdu=shavet;            
            xmp=xm;            
            nvd=nv0;
           end
%           'second', keyboard
           radii.a=[rapa(pu_p,sa); ailp; rapa(pu_d,sa)];
           radii.b=[rapb(pu_p,sa); bilp; rapb(pu_d,sa)];
           radii.c=[rapc(pu_p,sa); ailp*0; rapc(pu_d,sa)];
           shadd=ones(size(ailp))*shapar;
           shadd(1,:)=0;
           shadd(end,:)=0;           
           shavet=[shapdu(pu_p,sa); shadd; shapdu(pu_d,sa)];
           xm=[xmp(pu_p,san); -10*ones(size(nilp)); xmp(pu_d,san)];
           radii.array=[radii.array(pu_p1,:); radii.array(pudum,:); radii.array(pu_d,:)];
           nv0=[nvd(pu_p,san); nilp; nvd(pu_d,san)];

          else
           radii.a=[radii.a(pu_p,1); ailp; radii.a(pu_d,1)];
           radii.b=[radii.b(pu_p,1); bilp; radii.b(pu_d,1)];
           radii.c=[radii.c(pu_p,1); ailp*0; radii.c(pu_d,1)];
           shavet=[shavet(pu_p,1); ones(size(ailp))*0; shavet(pu_d,1)];
%          ' dopo  reass', keyboard
           xm=[xm(pu_p,1); -10*ones(size(ailp)); xm(pu_d,1)];
           radii.array=[radii.array(pu_p1,1); radii.array(pudum,1); radii.array(pu_d,1)];
           nva=nilp;
           nva(:,2)=0;
           nv0=[nv0(pu_p,1:2); nva; nv0(pu_d,1:2)];
          end 

%          ' dopo vero  part_ass ', keyboard
%         ' dopo vero  part_ass ', keyboard
%          ' dopo vero  part_ass ', keyboard

   
   
end   