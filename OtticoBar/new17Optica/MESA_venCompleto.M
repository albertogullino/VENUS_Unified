%'entro in MESA', keyboard

 if isfield(radii,'mesa')==1

  iMESA=zeros(size(shavet(:,1)));
%Isha_me=1|P40 Ra_me=15|P41 Dr_me=3|P42 StDmin=.1|P43 Include_me=0|P44 n_ext=1.58
%Mesa.Isha_me=3;
%Mesa.Drme=3;
%Mesa.Rame=15;
%Mesa.StDmin=3;
%Mesa.IMesa=1;
%Mesa.n_ext=1.67;



  Isham= pames{1};
  Ram  = pames{2};
  Drm  = pames{3};
  StDm = pames{4};
  Inclm= pames{5};
  n_ext= pames{6};
  ipar0= pames{7};
  ipme = pames{8};
  
%  Isham= pames{1};   forma pareti mesa: 0, verticali, 1: retta, n: polinomio grado n 
%  Ram  = pames{2};   raggio top (minimo)
%  Drm  = pames{3};    allargamento relativo base
%  StDm = pames{4};   n. discrettizzazioni
%  Inclm= pames{5};  non utilizzato
%  n_ext= pames{6};
%  ipar0= pames{7};
%  ipme = pames{8};   puntatore iniziale e finale zona mesa
  
  if ~exist('dilu')
   dilu=[];
  end
  if length(dilu)>0
  delme=length(dilu)-ipme(1)+1;
  else
  delme=0;
  end
  pume=(ipme(1):ipme(2))+delme;
  pume=pume-1;
%  'pume', keyboard
  if exist('Ip_grad')==1
  if length(find(pume==Ip_grad))==1
   pume=ipme(1):ipme(2)+Npadd_gra;
  end
  end
  iMESA(pume)=1;

  fipard=find(abs(fst(:,2))>1 & iMESA==1);
  difi=find(diff(fipard)>1);
  if length(difi)>0
   fipar=fipard(1:difi-1);
  else
   fipar=fipard;
  end



  raa=radii.a;
  rac=radii.c;
  naa=nv0;

Isham
%' Isham in mesa', keyboard

  if Isham==0   %pareti verticali
%'qui mesa verticale', keyboard
         for kr=pume
          rai=raa(kr,:);
          fao=find(rai>0);
          if length(fao)==0
           fao=0;
          end
          raa(kr,fao+1)=Ram;
          rac(kr,fao+1)=0;
          naa(kr,fao+2)=n_ext;
          shavet(kr,fao+1)=1;
         end
         radii.a=raa;
         radii.b=raa;
         radii.c=rac;
         nv0=naa;
%         for kp=pume
%          nau=(kp,:);
%          fiv=find(nau~=0);
%         end
         reassign
%'qui dopo mesa verticale', keyboard
     return
   end

  DTi=(abs(fst(:,2)).*dv.*iMESA);
  DT=sum(DTi)/1000;
  DTt=DT*1000;
  lox=Drm;
  DTip=(dv.*iMESA);
  DT_pair=sum(DTip(fipar))/1000;
  x=linspace(0,DT,100);
  
  Isham
%' Isham in mesa', keyboard
  
  if Isham>=1 %pareti retta
   pow=Isham;
   powi=1/Isham;
   Nd=fix(Drm/StDm);
   Del=DT;
   apar=(Del)/lox^pow;
   y=abs(x/apar).^powi;
   xd0=abs(linspace(0,max(abs(x)),Nd+1));
   xd=(xd0(1:end-1)+xd0(2:end))/2;
   yd=(xd/apar).^powi;
   xd=[xd0(1:end-1); xd0(2:end)];
   xd=reshape(xd,prod(size(xd)),1);
   yd=[yd; yd];
   yd=reshape(yd,prod(size(yd)),1);
   xt=[flipud(xd); ];
   yt=[flipud(yd); ];
   mabsy=max(abs(y));
   y=mabsy-y;
   yt=mabsy-yt;
   xd=max(x)-x;
   xdt=max(xt)-xt;
   if ifp==-10
    figure, plot(y,xd,yt,xdt,'.-')
    pausak
   end
   ipas_res=0;

    apart=apar*1000;
     idiscr_pair=1;
     xa=0;
     ya=0;
     na=[];
     ra=[];
     sa=[];
     yav=0;
     dvu=[];
     rvu=[];
     rep=[];
     rep1=[];
     ic=1;
     k=pume(ic);
     if ipver==1
      if ifp==-10
       hc=figure;
       set(hc,'pos',[452   271   560   420])
      end
     end
     while ic<=length(pume)
      k=pume(ic);
      its=fst(k,1);
      if its<=1
       ic=ic+1;
       xd=dv(k);
       nd=nv0(k,:);
       shd=shavet(k,:);
       rd=radii.a(k,:);
       xs=xa(end);
       xsf=DTt-xs;
       ys=Drm-abs((xsf)/apart).^powi;
       ys_sav=ys;
       if(abs(ys-ya(end)))<StDm/2
        if ipver==1
         'diff'
          abs(ys-ya(end))
          pausak
         end
        ys=ya(end);
       else
%        ' setto DU', pausak
       end
       yav=[yav; ys_sav];
       ya=[ya; ys];
       xa=[xa; xs+xd];
       dvu=[dvu; xd];
       rvu=[rvu; ys];
       rep=[rep; 1];
       rep1=[rep1; 0];
       na=[na; nd];
       ra=[ra; rd];
       sa=[sa; shd];
       if ipver==1
         xd
         xa
         ys_sav
         '[xs ys]=',[xs ys]
        if ifp==-10
         figure(hc), plot(xa,ya,'r.')
         ' its =1 ', pausak
        end
       end
      else
       kv=k:k+its-1;
       ic=ic+its;
       xd=dv(kv);
       nd=nv0(kv,:);
       shd=shavet(kv,:);
       rd=radii.a(kv,:);
       repi=1;
%       repis=repi;
       for kkm=1:fst(k,2)
         ireset=1;
        xs=xa(end);
        xs=xs+cumsum(xd);
        if ipver==1
         kkm
         xd
         xa
         ys_sav
        end
        xa=[xa; xs];
        xsf=DTt-xs;
        ys=Drm-abs((xsf)/apart).^powi;
        ys_sav=ys;
        yav=[yav; ys_sav];
        if (max(abs(ys-ya(end))))<StDm/2
         if ipver==1
          'diff'
          abs(ys-ya(end))
          pausak
         end
         ireset=0;
         ys=ya(end)*ones(size(ys));
         repi=repi+1;
%         repis=repi;
        else
         if StDm/DT_pair>0.7
%          ys=yav(end)*ones(size(ys))+StDm;
          ys=yav(end)*ones(size(ys));
         end
%        ' setto DU mirro', pausak
        end
        ya=[ya; ys];
        if ipver==1
         [xs ys]
         if ifp==-10
          figure(hc), plot(xa,ya,'r.')
%         '[repi repis]'
%         [repi repis]
         'repi=',repi

          pausak
         end
        end
        if (ipas_res==0 & kkm==fst(k,2) )
         na=[na; nd];
         ra=[ra; rd];
         sa=[sa; shd];
         dvu=[dvu; xd];
         rvu=[rvu; ys];
         rep=[rep; fst(kv,2)];
         rep1=[rep1; fst(kv,1)];
        end
        if ireset==1
         ipas_res=1;
         na=[na; nd];
         ra=[ra; rd];
         sa=[sa; shd];
         dvu=[dvu; xd];
         rvu=[rvu; ys];
         rep=[rep; ones(size(ys))*repi];
         rep1=[rep1; ones(size(ys))*fst(k,1)];
         repi=1;
        end
       end
%        if repi>1
%         dvu=[dvu; xd];
%         rvu=[rvu; ys];
%         rep=[rep; ones(size(ys))*repi];
%        end
        if ipver==1
         if ifp==-10
         figure(hc), plot(xa,ya,'r.')
         ' its >1 ', pausak
         end
        end
      end
     end

     rmes=Ram+rvu;
     rept=[rep1 rep];
     nmes=rmes;
     nmes(:)=n_ext;
     shadd=sa;

   xad=xa;
   if ifp==-10
   figure, plot(ya,xa,'.-',yav,xa,'r.-'), title(' discretizzato - vero')
   pausak
   ' ferma ', pausak
   end
   
    rmeto=ra;
  shmeto=shadd;
  sha0=1;
%  ' metao', keyboard
  for k=1:length(ra)
   rl=rmeto(k,:);
   fi0=find(rl==0);
   if length(fi0)==0
    rmeto(k,end+1)=rmes(k);
    shmeto(k,end+1)=sha0;
   else
    rmeto(k,fi0(1))=rmes(k);
    shmeto(k,fi0(1))=sha0;
   end
%   [k shmeto(k,:)]
%   pausak
  end


  nilu=na;
  for k=1:length(na)
   rl=nilu(k,:);
   fi0=find(rl==0);
   if length(fi0)==0
    nilu(k,end+1)=n_ext;
   else
    nilu(k,fi0(1))=n_ext;
   end
  end


          dilu=dvu;
          ailu=rmeto;
          ipudv=pume;

          pu_p=1:ipudv(1)-1;
          pu_p1=1:ipudv(1);
          pu_d=ipudv(end)+1:length(dv);
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

          dv=[dv(pu_p); dilu; dv(pu_d)];
          ipar=[ipar(pu_p1,:,:); dua; ipar(pu_d,:,:)];
          anyf=[anyf(pu_p); duze; anyf(pu_d)];
          fst=[fst(pu_p,:); rept; fst(pu_d,:)];

          iau_add=zeros(length(dilu)-1,size(iauto,2));
          iau_add=[iau_add; iauto(pu_d(1)-1,:)];

 %         iauto=[iauto(pu_p,:); iau_add; iauto(pu_d,:) ];          
%          iauto=[iauto(pu_p,:); [zeros(size(dilu)) zeros(size(dilu))]; iauto(pu_d,:) ];
          ifield=[ifield(pu_p,:); zeros(size(dilu)); ifield(pu_d,:)];
          pudum=ones(size(dilu(1:end-1)))*(pu_p1(end)+1);
          sail=size(ailu);
          sa=1:sail(2);
          sar=size(radii.a);
          if sa(end)>sar(2)
           radii.a=[radii.a zeros(length(radii.a),sa(end)-sar(2))];
           radii.b=[radii.b zeros(length(radii.a),sa(end)-sar(2))];
           radii.c=[radii.c zeros(length(radii.a),sa(end)-sar(2))];
           shavet=[shavet zeros(length(shavet),sa(end)-sar(2))];
           xm=[xm zeros(length(radii.a),sa(end)-sar(2))];
           nv0=[nv0 zeros(length(radii.a),sa(end)-sar(2))];
          end
          said=size(nilu);
          san=1:said(2);
          radii.a=[radii.a(pu_p,sa); ailu; radii.a(pu_d,sa)];
          radii.b=[radii.b(pu_p,sa); ailu*0; radii.b(pu_d,sa)];
          radii.c=[radii.c(pu_p,sa); ailu*0; radii.c(pu_d,sa)];
          shavet=[shavet(pu_p,sa); shmeto; shavet(pu_d,sa)];
          xm=[xm(pu_p,san); -10*ones(size(nilu)); xm(pu_d,san)];
          radii.array=[radii.array(pu_p1,:); radii.array(pudum,:); ...
                       radii.array(pu_d,:)];
          nv0=[nv0(pu_p,san); nilu; nv0(pu_d,san)];
          clear ipudv

          %fishav=find(shavet<0 & abs(shavet)<=4);
          %iGRAD=zeros(size(shavet));
          %if length(fishav)>0
          % iGRAD(fishav)=1;
          %end 
   
  end


  
 if Isham==-1 %pareti X
  
  %' qui Isham Xcsel', keyboard
   unpack_x
   fiox=find(ntoe==1.6);
   dvmul=[dto];

   eth_angle=Ps.Mesa.eth_angle; 
   Rmaxt=Ps.Mesa.Ratop;
   Ram=Ps.Mesa.Ramin;   
   zox=sum((dvmul(1:fiox-1)));
   zoxi=sum((dvmul(fiox:end)));
   dvtot=cumsum(flipud(dvmul));
   dvtot1=cumsum(flipud(dvmul))-zoxi;
   ratot=Ram+abs(dvtot1)/1000*tan((90-eth_angle)/180*pi);
   Rdis=Ps.Mesa.Rdis;
   Lnomet=Ps.Mesa.Lnomet;
   znomet=zoxi-Lnomet*cos((90-eth_angle)/180*pi);
   Next=Ps.Mesa.next;
   Nextc=Ps.Mesa.nextc; 
   if isfield(Ps.Mesa,'dext')==1
    Dext=Ps.Mesa.dext/1000;
   else
    Dext=0;
   end
   Rdis_start=Ps.Mesa.Rdis_start;
   ra_val=[Ram+Rdis_start:Rdis:max(ratot)];
   ra_dis=ratot;
   zmetal=zeros(size(ratot));
   for krag=2:length(ra_val)
    fi=find(ratot<ra_val(krag) & ratot>=ra_val(krag-1));
    ra_dis(fi,1)=ra_val(krag-1);
    if ra_val(krag-1)>Rmaxt
     fize=find(dvtot(fi)>zoxi);
     if length(fize)>0
      fizev=fi(fize);
      ra_dis(fizev,1)=Rmaxt;     
      zmetal(fizev,1)=1;     
     end
    end 
   end
   ra_dis(end)=Rmaxt;
   fizet=find(dvtot<znomet);
   zmetal(fizet,1)=1;
   zmetal(end,1)=1;
   fibot=find(ra_dis>Rmaxt+1);
%   fibot=find(ra_dis>Ram+1);
%   fibot=find(ra_dis>Ram+.5);
%   fibot=find(ra_dis>Ram+.2);
%   fibot=find(ra_dis>Ram+.3);
%   fibot=find(ra_dis>Ram+.4);
   zmetal(fibot,1)=-1;
   flzmet=flipud(zmetal);
   fimet=find(zmetal==1);
   

%   if ifp==-10
    figure, plot(ratot,(dvtot/1000),'.',ra_dis,(dvtot/1000),'r.',ra_dis(fimet),dvtot(fimet)/1000,'c.'),
    axis([0 max(ratot) 0 10])
    axis equal   
    pausak
%   end 

   


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


   pma=max(fsto);
   fima1=find(fsto<pma & fsto>1);
   pma1=fsto(fima1(1));
   fi_su=find(fst(:,2)==pma);
   paper=fst(fi_su(1),1);
   fi_in=find(fst(:,2)==pma1);
   paper1=fst(fi_in(1),1);   
   
   fiparu=find(flzmet==1 & fsto==pma);
   nPsup=floor(length(fiparu)/paper);
   fiparu=find(flzmet==-1 & fsto==pma1);
   nPinf=floor(length(fiparu)/paper1);
   
      fi_sut=find(fsto==pma);
      fi_int=find(fsto==pma1);
   per_sup=fi_sut(1)+[0:paper-1];
   per_inf=fi_int(end)-[paper1-1:-1:0];
   
   fsto1=ones(size(fsto));
   fsto1(per_sup)=paper;
   fsto(per_sup)=nPsup;
   for k=1:nPsup-1
      fsto(per_sup+k*paper)=0;
   end
   fsto1(per_inf)=paper1;   
   fsto(per_inf)=nPinf;   
   for k=1:nPinf-1
      fsto(per_inf-k*paper)=0;
   end

      fi_sut=find(fsto==pma);
      fsto(fi_sut)=1;
      fi_int=find(fsto==pma1);   
      fsto(fi_int)=1;

fivalidi=find(fsto~=0);

'qui', pausak


          flzmet=flzmet(fivalidi);
          dv=dvmul(fivalidi);

          nto=nto(fivalidi);
          rto=rto(fivalidi);
          ntoe=ntoe(fivalidi);
%          ipar=[ipar(pu_p1,:,:); dua; ipar(pu_d,:,:)];
          anyf=zeros(size(dv));
          ifield=zeros(size(dv));
          fst=[fsto1(fivalidi) fsto(fivalidi)];
          ratot=ra_dis;
          ravx=[flipud(ratot)]; 
          ravx=ravx(fivalidi,:);
          raddi=[];
          if Dext>0
           for irag=1:length(Dext)                        
            raddi=[raddi (ravx+sum(Dext(1:irag))).*flzmet];
           end 
          end
           ravx=[ravx raddi];

          
          nvx=[ones(size(nto))*Next];
          
          for kla=1:length(ravx)
           fi=find(rto(kla)>0);
           rvec=[];
           nvec=[nto(kla)];
           if length(fi)>0
            rvec=rto(kla); 
            nvec=[ nvec ntoe(kla)]; 
           end
           if flzmet(kla)==1
            rtudu=[rvec ravx(kla,:)];
            ntudu=[nvec nvx(kla,:)];
           elseif flzmet(kla)==0
            rtudu=[rvec ravx(kla,1)];
            ntudu=[nvec Nextc];  
           else 
            rtudu=rvec;
            ntudu=nvec; 
           end
           rtu(kla,1:length(rtudu))=rtudu;

           ntu(kla,1:length(ntudu))=ntudu;  
           radii.array{kla}{1}=0;
          end 
          raver=rtu;
          nver=ntu;
          
           radii.a=raver;
           radii.b=raver*0;
           radii.c=raver*0;
           shavet=[ones(size(raver))];
           xm=cto(fivalidi);
           nv0=nver;
           iauto=iauo(fivalidi,:);
           iXcsel=1;
           
  end




    reassign
     dvp=[dv(1:end-1) dv(2:end)];
     dvp=reshape(dvp,prod(size(dvp)),1);




end


'fine in MESA', keyboard