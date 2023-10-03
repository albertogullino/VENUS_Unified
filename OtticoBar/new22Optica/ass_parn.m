%' ass_parn', keyboard
ipas_grat=0;
ipver=1;   % mostra tutto per verifica

%' in ass_parn ', keyboard

     fish=find(shavet==6 | shavet==-4);
   
     if length(fish)>0  
       par=radii.array{fish};
%     'AssParn', keyboard
       ngrating=nv0(fish,:);
       ngrating1=ngrating(1);
       ngrating2=ngrating(2);
       if  length(par)==13
        if length(par{12})>0
         ngrating1=par{13};
         par_grat.th=par{12};
        end
       end
       ngrating2=ngrating(2)*ones(size(ngrating1));
       par_grat.r_in=nv0(fish-1,1);
       par_grat.r_out=nv0(fish+1,1);
       par_grat.n1=ngrating1;
       par_grat.n2=ngrating2;
       par_grat.r1=ngrating1;
       par_grat.r2=ngrating2;       
       
     
       
       n1g=ngrating1;
       n2g=ngrating2;
       period=par{5};
       t1=par{6};
       t2=period-t1;
       DC=t1/period;
       par_grat.itetm=3;
       par_grat.px=period;
       par_grat.per=period;
       par_grat.DC=DC;
       if isfield(Pf,'NModi')==1
        NModi=Pf.NModi;
       else
        NModi=11;
       end
       par_grat.NModi=NModi;
       itetm=3;
       icarico=0;   % calcola T reticolo caricato sui G_i
       dret=L_i(fish);
       %' cont grat', keyboard
       ibast=fish;
     else 
       ibast=[];
       dret=0;
     end  
 

%'MESA', keyboard

 iset_MESA=0;
 iset_GRAD=0;
 if isfield(radii,'mesa')==1
  if isfield(Ps,'Mesa')==1
   if isfield(Ps.Mesa,'Isha_me')==1
    if Ps.Mesa.Isha_me>=0
     iset_MESA=1;
     %PAME
     pames=radii.mesa;
      if isfield(Ps.Mesa,'Isha_me')
       pames{1}=Ps.Mesa.Isha_me;
      end 
      if isfield(Ps.Mesa,'Rame')
       pames{2}=Ps.Mesa.Rame;      
      end
      if isfield(Ps.Mesa,'Drme')
       pames{3}=Ps.Mesa.Drme;      
      end   
      if isfield(Ps.Mesa,'StDmin')
       pames{4}=Ps.Mesa.StDmin;      
      end         
      if isfield(Ps.Mesa,'n_ext')
       pames{6}=Ps.Mesa.n_ext;          
      end 
       
%  Isham= pames{1};   forma pareti mesa: 0, verticali, 1: retta, n: polinomio grado n 
%  Ram  = pames{2};   raggio top (minimo)
%  Drm  = pames{3};    allargamento relativo base
%  StDm = pames{4};   n. discrettizzazioni
%  Inclm= pames{5};  non utilizzato
%  n_ext= pames{6};
%  ipar0= pames{7};
%  ipme = pames{8};   puntatore iniziale e finale zona mesa

    end
   end
  end
 end


    

%'entro asspaer MESA', keyboard

 s=size(par_in);
 lpar=s(1);
 iset_DOE=0;
 shavet_le=shavet(:,1);
 fishav=find(shavet<0 & abs(shavet)<=4 & shavet~=-4);
 if length(fishav)>1
  ' errore, errore in ass_par: non piu di uno strato graded '
  ' programma valido solo per uno strato graded '
  keyboard
 end
 if length(fishav)==1
  if isfield(Ps,'Grad')==1
   fna=fieldnames(Ps.Grad);
   for ifna=1:length(fna)
    pav{ifna}=eval(['Ps.Grad.',fna{ifna}]);
   end
  else
   pav=radii.array{fishav};
  end
  Ishag=pav{1};
%  ' iasi',
%  keyboard
  if Ishag>0
   iGRAD=zeros(size(shavet));
   iset_GRAD=1;
   iGRAD(fishav)=1;
  else
   iset_GRAD=0;
  end
 end


 if iset_GRAD==1
  pav=radii.array{fishav};
  PAlgrad.Isha=pav{1};
  PAlgrad.Dr=pav{2};
  PAlgrad.Nr=pav{3};
  PAlgrad.Iapp=pav{4};

  dus1=reshape(ipar(fishav,1,:),9,1);
  dus2=reshape(ipar(fishav,2,:),9,1);
  fid=find(dus2==-2);
  if length(fid)==1
   PAlgrad.Isha=par_in{dus1(fid)};
  end
  fid=find(dus2==-3);
  if length(fid)==1
   PAlgrad.Dr=par_in{dus1(fid)};
  end
  fid=find(dus2==-4);
  if length(fid)==1
   PAlgrad.Nr=par_in{dus1(fid)};
  end
  fid=find(dus2==-5);
  if length(fid)==1
   PAlgrad.Iapp=par_in{dus1(fid)};
  end
 end

 idoe=[0 0];

 if lpar>0

  sip=size(ipar);
  riga=[1:sip(1)]';
  righe=riga*ones(1,sip(3));
  ip1=reshape(ipar(:,1,:),sip(1)*sip(3),1);   % par. #
  ip2=reshape(ipar(:,2,:),sip(1)*sip(3),1);   % par. type
  ip3=reshape(ipar(:,3,:),sip(1)*sip(3),1);   % trans. sect.
  ip4=reshape(righe,sip(1)*sip(3),1);         % layer
  fi1=find(ip1~=0);
  [dup,fi1p]=sort(ip1(fi1));
  fi1=fi1(fi1p);

  fiu=fi1;

  npar=length(find(diff([0; dup])~=0));

% Bisogna fare il MESA dopo

lfor=length(par_in(:,1));
ipar_all=1:lfor;
ipar_altri=1:lfor;

if iset_MESA==1
  iparme=radii.mesa{7};
  sipa=size(iparme);
  if length(sipa)==3
   iparmea=reshape(iparme(1,1,:),1,sipa(3));
   iparmea2=abs(reshape(iparme(1,2,:),1,sipa(3)))-1;
  else
   iparmea=iparme(1,1,1);
   iparmea2=abs(iparme(1,2,1))-1;  
  end
  iparal=ipar_altri;
  for ii=1:length(iparmea)
   fid=find(iparmea(ii)==iparal);
   iparal(fid)=-1;
  end
  fid=find(iparal~=-1);
%  ipar_altri=ipar_altri(fid);
else
 iparmea=[];
end



% Bisogna fare il MESA dopo
% Faccio il resto (NO MESA)
%  if lpar<npar
%  if lpar>npar
  if lpar==-10

   disp(' error: par_in does not match the number of parameters ')
   keyboard
   return
  else

  lfor=length(par_in(:,1));

%   for ip=ipar_altri %tutti in realta'
   for ip=ipar_all
%    ip, pausak
     if length(find(iparmea==ip))==1
%      ' setto valore mesa', pausak

%      switch par
%       case 'Isha_me='
%       nupar=-2;
%       case 'Dr_me='
%       nupar=-3;
%       case 'StDmin='
%       nupar=-4;
%       case 'Include_me='
%       nupar=-5;
%       case 'Ra_me='
%       nupar=-6;
%       case 'n_ext='
%       nupar=-7;
%      end
%  Isham= pames{1};
%  Ram  = pames{2};
%  Drm  = pames{3};
%  StDm = pames{4};
%  Inclm= pames{5};
%  n_ext= pames{6};
%  ipar0= pames{7};
%  ipme = pames{8};      
      
      pame=pames{7}(2);
      nuparmes=[-2 -3 -4 -5 -6 -7];
      fime=find(nuparmes==pame);
      nuordmes=[ 1  3  4   5  2  6];
      pu_me=nuordmes(fime);
%      pu_me=ip+1-iparmea(1);
      PAram=par_in{ip,1};
      pames{pu_me}=PAram;
%      'fine setto valore mesa', pausak

     end

      fip=find(ip1(fiu)==ip);
%      pausak

      lpai=length(par_in{ip,1});

      if length(fip)>0 & lpai>0
       pup=fiu(fip);
       pur=ip4(pup);
       puc=ip3(pup);
       ipad=ip2(pup);
       ipa1=ip1(pup);
       ipav=ipad;
%'       [ipa1 ipad puc pur]'
%       [ipa1 ipad puc pur]
%       pausak

%    ipa=1:  thickness
%    ipa=2:  y-axis
%    ipa=3:  x-axis
%    ipa=4:  shape param (when rhombus)
%    ipa=5:  ref. index
%    ipa=6:  array centers
%    ipa=7:  type of shape
%    ipa=8:  # of repetitions
%    ipa=9:  additional circle radius (when grating)
%    ipa=10:  additional Rext circle radius (when grating)

%     if ip>5
%      disp('ass_par, ipav'), ipav, keyboard
%      pausak
%     end
icpa=0;
       for ipa=ipav'
       icpa=icpa+1;
%       '[ip, ipa]',         [ip, ipa], pausak

        if ipa>0

         if ipa==1
%          dv(pur)=par_inn(ip);
          dv(pur)=par_in{ip,1};
         elseif ipa==2
%          radii.a(pur,puc)=par_in(ip,1);
          for isr=1:length(pur)
           radii.a(pur(isr),puc(isr))=par_in{ip,1};
          end
         elseif ipa==3
%          radii.b(pur,puc)=par_in{ip,1};
          for isr=1:length(pur)
           radii.b(pur(isr),puc(isr))=par_in{ip,1};
          end
         elseif ipa==4
%          radii.c(pur,puc)=par_in{ip,1};
          for isr=1:length(pur)
           radii.c(pur(isr),puc(isr))=par_in{ip,1};
          end
         elseif ipa==5
%          nv0(pur,puc)=par_in{ip,1};
          for isr=1:length(pur)
           if ipad(isr)>0
            nv0(pur(isr),puc(isr))=par_in{ip,1};
           end 
          end
%          'ipa = 5, ref Index', keyboard
         elseif ipa==10
%          nv0(pur,puc)=par_in{ip,1};
          %'ipa = 10, ref Index', keyboard          
          
          for isr=1:length(pur)
           nv0(pur(isr),puc(isr))=real(nv0(pur(isr),puc(isr)))-j*par_in{ip,1}*50/k0;
          end
%          'ipa = 10, ref Index', keyboard          
         elseif ipa==6
          'ipa = 6', keyboard
          nce=length(par_in{ip,3});
          radii.array{pur,2}=nce;
          for kdu=3:7
           radii.array{pur,kdu}=par_in{ip,kdu};
          end

%         radii.array{pcdu,2}=ipcen-1;
%         radii.array{pcdu,3}=cemem;
%         radii.array{pcdu,4}=shapes;
%         radii.array{pcdu,5}=Ry;
%         radii.array{pcdu,6}=Rx;
%         radii.array{pcdu,7}=Delta;

         elseif ipa==7
%          shavet(pur)=par_inn(ip,1);
%          'shavet', keyboard
%          shavet(pur)=sign(shavet(pur))*shaf(par_in{ip,1});
          shame=shavet(pur(1));
%          ' keyboard', keyboard
          shavet(pur)=shaf(par_in{ip,1});
          fifi=find(shavet==shame);
          if length(fifi)>0
           shavet(fifi)=shaf(par_in{ip,1});
          end
          if par_in{ip,1}==0
           radii.a(pur,puc)=0;
          end
         elseif ipa==8
%          fst(pur,2)=par_inn(ip,1);
          fst(pur,2)=par_in{ip,1};
         end

        else  %ipa<0


        if iset_GRAD==1
         if iGRAD(pur)==1
           PAram=par_in{ip,1};
          if ipa==-2
           PAlgrad.Isha=PAram;
          elseif ipa==-3
           PAlgrad.Dr=PAram;
          elseif ipa==-4
           PAlgrad.Nr=PAram;
          elseif ipa==-5
           PAlgrad.Iapp=PAram;
          end
         end
        end
        if ifp==-10
%         ' < 0 ', keyboard
        end
         shapuv=abs(shavet(pur,1));
         purv=pur;
%         pur
%         shapuv
%         ' < 0 ', keyboard

         for kassegna=1:length(purv)
          purloc=purv(kassegna);
          shpu=shapuv(kassegna);
%          shpu, pausak
           if shpu<6
            ik=0;
            for kpu=purloc'
             ik=ik+1;
             radii.array{kpu,puc(ik)}{abs(ipa)}=par_in{ip};
            end
           else
            kpu=purloc;
            if shpu==7
%            ip, kpu
%             ' PAle ', keyboard
             PAle{puc}=radii.array{kpu}{11};
             PAth{puc}=radii.array{kpu}{12};
             PAnr{puc}=radii.array{kpu}{13};
%             ' ipudv ', keyboard
             if ~exist('ipudv')
%              ipudvf=find(radii.a(kpu:end,1)>0);
              ipudvf=find(shavet_le==7);
%              'ipudvf', keyboard
%              'ipudvf', keyboard
%              'ipudvf', keyboard
%              ipudv=kpu-1+ipudvf;
              ipudv=ipudvf;
             end
%             [ipa ip]
%            ' aass', pausak
             if ipa==-3
              PAth{puc}=par_in{ip};
             elseif ipa==-4
              PAnr{puc}=par_in{ip};
             elseif ipa==-5
              PAle{puc}.H=par_in{ip};
%              'ver H', keyboard
             elseif ipa==-6
              PAle{puc}.D=par_in{ip};
             elseif ipa==-7
              PAle{puc}.Ndis=par_in{ip};
             elseif ipa==-8
              PAle{puc}.NlAR=par_in{ip};
             elseif ipa==-10
              PAle{puc}.Rel=par_in{ip};
             elseif ipa==-11
              PAle{puc}.Npair=par_in{ip};
             elseif ipa==-12
              PAle{puc}.Rflat=par_in{ip};
             elseif ipa==-13
              PAle{puc}.Nrel=par_in{ip};
             elseif ipa==-14
              PAle{puc}.Rm_rel=par_in{ip};
             elseif ipa==-15
              PAle{puc}.Rax=par_in{ip};
             elseif ipa==-16
              PAle{puc}.Mis=par_in{ip};              
             elseif ipa==-17
              PAle{puc}.pos=par_in{ip};     
             elseif ipa==-18
              PAle{puc}.thg=par_in{ip};     
             elseif ipa==-19
              PAle{puc}.thga=par_in{ip};                   
             elseif ipa==-20
              PAle{puc}.LA=par_in{ip};                   
             elseif ipa==-21
              PAle{puc}.d=par_in{ip};    
             elseif ipa==-22
              PAle{puc}.orien=par_in{ip};                  
             elseif ipa==-23
              PAle{puc}.Base=par_in{ip};  
             elseif ipa==-24
              PAle{puc}.Alt=par_in{ip};   
             elseif ipa==-25
              PAle{puc}.Sth=par_in{ip};                 
%              ' Sth', keyboard
             end
             radii.array{kpu}{11}=PAle{puc};
             radii.array{kpu}{12}=PAth{puc};
             radii.array{kpu}{13}=PAnr{puc};
%'lente', keyboard

            elseif shpu==6 
            
              if isfield(Ps,'igraef_new')==1
               igraef_new=Ps.igraef_new;
              else 
               igraef_new=-1;
              end
              pa_temp=radii.array{kpu}{11};
%              if igraef_new==2
            if isstruct(pa_temp)==0
               puci=puc(kassegna);
%            ' aass', pausak
%'vedo ipa in HCG', ipa, pausak
             if ipa==-3
              par_grat.th(puci)=par_in{ip};
             elseif ipa==-5
              par_grat.px=par_in{ip};
              par_grat.per=par_in{ip};
             elseif ipa==-6
              par_grat.DC=par_in{ip};
             elseif ipa==-4
%              ' par_gr', keyboard
              par_grat.r1(puci)=par_in{ip}+imag(par_grat.r1(puci));
              par_grat.n1(puci)=par_in{ip}+imag(par_grat.r1(puci));

             elseif ipa==-40
              par_grat.r1(puci)=-j*1e4/(2*k0)*par_in{ip}+real(par_grat.r1(puci));
              par_grat.n1(puci)=-j*1e4/(2*k0)*par_in{ip}+real(par_grat.r1(puci));
             end 
%'dopo ipa in HCG', ipa, pausak             
%	               ngrating1=par{13};
%	               par_grat.th=par{12};

%
             
            else
%             par_in{ip}
%             keyboard
                         if ipa==-30
	                  pa_temp.th(puc)=par_in{ip};
	                  dv(pur)=sum(pa_temp.th);
	                 elseif ipa==-41
	                  pa_temp.dcv(puc)=par_in{ip};
	                 elseif ipa==-5
	                  pa_temp.Pe=par_in{ip};
                         end 
                         radii.array{kpu}{11}=pa_temp;
%               ' HCG par_grat ', keyboard         

            end %is Struct

            elseif shpu==8
% in fipar.m dentro lay_tran
%       case 'RoC='
%       nupar=-5;
%       case 'R='
%       nupar=-6;
%       case 'DC='
%       nupar=-7;
%       case 'Pe='
%       nupar=-8;   
%       case 'Ns='
%       nupar=-9;  

             if ipa==-5
              iset_DOE=1;
              pe=par_in{ip};
              idoe(1)=1;
             elseif ipa==-6
              iset_DOE=1;
              D=par_in{ip};
              radii.array{kpu}{12}=D;              
              radii.a(kpu)=D;              
              idoe(2)=1;
             elseif ipa==-7
              DC=par_in{ip};
              PG=radii.array{kpu}{13};
              PG.DC=DC;
              radii.array{kpu}{13}=PG;              
              radii.array{kpu}{7}=DC;
             elseif ipa==-8
              PE=par_in{ip};
              PG=radii.array{kpu}{13};
              PG.period=PE;
              radii.array{kpu}{13}=PG;              
             end
             
           %'DOE', keyboard
            elseif shpu==9 
               radii.array{kpu}{abs(ipa)}=par_in{ip};
               PApa=radii.array{kpu};             
%             'ipa, ip', [ipa ip]
%             'son chgi', keyboard
            end
           end
%           [shpu ip]
%           ' qui asspar ', pausak
         end %shapu
        end

       end  %for ipa
%     end  % mesa
%     disp(' loo 2 '), pausak

    end
   end
   if exist('pa_temp')
           if isstruct(pa_temp)==1
%       'HCG', keyboard
       'HCG multiplo'
             par_grat.DC=pa_temp.dcv(1);
             par_grat.n1=pa_temp.nr1;
             par_grat.r1=pa_temp.nr1;
             par_grat.n2=pa_temp.nr2;
             par_grat.r2=pa_temp.nr2;       
             par_grat.px=pa_temp.Pe;
             par_grat.per=pa_temp.Pe;
             par_grat.th=pa_temp.th;
            end 
      else      
       if fish>0
 	kpu=fish;
             radii.array{kpu}{5}=par_grat.per;
             radii.array{kpu}{6}=par_grat.per*par_grat.DC;
             if isfield(par_grat,'th')==1
              radii.array{kpu}{12}=par_grat.th;
              radii.array{kpu}{13}=par_grat.n1;
             end
	end
       end	 
%   ' fine parziale parametri ', keyboard
 %  ' fine parziale parametri ', keyboard
 %  ' fine parziale parametri ', keyboard

      if iset_GRAD==1
       lox=dv(fishav);
       Ndis=PAlgrad.Nr;
       Del= PAlgrad.Dr;
       Isha=PAlgrad.Isha; %2 parab, 1 retta, 0 cost
       if Isha==2
        Nd=ceil(Ndis/2);
        x=linspace(-lox/2,lox/2,50);
        apar=(4*Del)/lox^2;
        y=apar*x.^2;
        xd0=linspace(0,lox/2,Nd+1);
        xd=(xd0(1:end-1)+xd0(2:end))/2;
        yd=apar*xd.^2;
        xd=[xd0(1:end-1); xd0(2:end)];
        xd=reshape(xd,prod(size(xd)),1);
        yd=[yd; yd];
        yd=reshape(yd,prod(size(yd)),1);

        xt=[-flipud(xd(2:end)); xd(2:end)];
        yt=[flipud(yd(2:end)); yd(2:end)];
        if ifp==-10
         figure, plot(x,y,xt,yt,'.-')
         pausak
        end
        ydu=diff(xt);
        fdu=find(ydu>0);
        y0=radii.a(fishav);
        xt_lay=ydu(fdu);
        yt_lay=yt(fdu)+y0;
       elseif Isha==1
        Nd=Ndis;
        x=linspace(0,lox,50);
        apar=-(Del)/lox;
        y=apar*x;
        xd0=linspace(0,lox,Nd+1);
        xd=(xd0(1:end-1)+xd0(2:end))/2;
        yd=apar*xd;
        xd=[xd0(1:end-1); xd0(2:end)];
        xt=reshape(xd,prod(size(xd)),1);
        yd=[yd; yd];
        yt=reshape(yd,prod(size(yd)),1);
        if ifp==-10
         figure, plot(x,y,xt,yt,'.-')
         pausak
        end
        ydu=diff(xt);
        fdu=find(ydu>0);
        y0=radii.a(fishav);
        xt_lay=ydu(fdu);
        yt_lay=yt(fdu)+y0;

       elseif Isha==0
       end

       if Isha>0

        if PAlgrad.Iapp==0
          dilu=xt_lay;
          Npadd_gra=length(dilu)-1;
          Ip_grad=fishav;
          ailu=yt_lay;
          nilu=ones(size(ailu))*nv0(fishav,:);
          ipudv=fishav;
%          itrle=find(diff(ipudvf)>1)
%          keyboard
%          if length(itrle)==0
%           itrle=length(ipudv);
%          else
%           shavet_le(ipudvf(1:itrle))=0;
%          end

          pu_p=1:ipudv(1)-1;
          pu_p1=1:ipudv(1);
%          pu_d=ipudv(itrle)+1:length(dv);
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
          shaox=shavet(fishav);
          n_sa=nv0;
          xm_sa=xm;
%'QUI diliu', keyboard

          dv=[dv(pu_p); dilu; dv(pu_d)];
          ipar=[ipar(pu_p1,:,:); dua; ipar(pu_d,:,:)];
          anyf=[anyf(pu_p); duze; anyf(pu_d)];
          fst=[fst(pu_p,:); [zeros(size(dilu)) ones(size(dilu))]; fst(pu_d,:)];
          iauto=[iauto(pu_p,:); [zeros(size(dilu)) zeros(size(dilu))]; iauto(pu_d,:) ];
          ifield=[ifield(pu_p,:); zeros(size(dilu)); ifield(pu_d,:)];
          pudum=ones(size(dilu(1:end-1)))*(pu_p1(end)+1);
          sail=size(ailu);
          sa=1:sail(2);
          said=size(nilu);
          san=1:said(2);
          radii.a=[radii.a(pu_p,sa); ailu; radii.a(pu_d,sa)];
          radii.b=[radii.b(pu_p,sa); ailu*0; radii.b(pu_d,sa)];
          radii.c=[radii.c(pu_p,sa); ailu*0; radii.c(pu_d,sa)];

          shavet=[shavet(pu_p,sa); ones(size(ailu))*shaox; shavet(pu_d,sa)];
          xm=[xm(pu_p,san); -10*ones(size(nilu)); xm(pu_d,san)];
          radii.array=[radii.array(pu_p1,:); radii.array(pudum,:); radii.array(pu_d,:)];
          nv0=[nv0(pu_p,san); nilu; nv0(pu_d,san)];
          clear ipudv
        else
          ailu_mean=yt_lay;
          dilu_mean=xt_lay;
          iKmean=1;
          radii.a(fishav)=-radii.a(fishav);
        end
       end  %Isha

      end

      if iset_DOE==1
         if prod(idoe(1:2))==1
          Dd=pe-D;
          for kdu=kpu'
           radii.array{kdu}{5}=D;
           radii.array{kdu}{6}=Dd;
          end
          iDOEver=1;
         else
          iDOEver=0;
          for kdu=kpu'
           radii.array{kdu}{12}=D;
           radii.a(kdu)=D;
          end 
%          ' errore ass_par 123 ', keyboard
         end

        if iDOEver==1  
         pe0=[D D+Dd];
         duth=pe0;
         k=1;
         while duth(end)<=Ry
          duth=[duth pe0+k*pe];
          k=k+1;
         end

         dup=[[0; duth(1:end-1)'] [duth(1:end-1)'; duth(end)+10]];
         val=zeros(size(dup));
         val(1:2:end,:)=1;
         ru=reshape(dup',1,prod(size(dup)));
         pu=reshape(val',1,prod(size(dup)));
         if ifp~=-4
         figure, plot(ru,pu,'.-r'), axis([0 ru(end) 0 1.1]); pausak
         end

         for kdu=kpu'
          radii.array{kdu}{10}=duth(1:end-1);
         end
        else
%' ver', keyboard
        end
       end  %iset_DOE

       if exist('ipudv')    %lente
        for klen=1:length(PAle)
          duth=PAth{klen};
          dunr=PAnr{klen};
          ha=PAle{klen}.H;
          ra=PAle{klen}.D;
          Ndisc=real(PAle{klen}.Ndis);
          UD=PAle{klen}.UD;
          Rel=PAle{klen}.Rel;
          Nrel=PAle{klen}.Nrel;
          if isfield(PAle{klen},'Npair')==1
           Npair=PAle{klen}.Npair;
          else
           Npair=0;
          end
          if isfield(PAle{klen},'Rflat')==1
           Rflat=PAle{klen}.Rflat;
          else
           Rflat=0;
          end
          if isfield(PAle{klen},'Rm_rel')==1
           Rm_rel=PAle{klen}.Rm_rel;
          else
           Rm_rel=0;
          end
          if isfield(PAle{klen},'Rint_rel')==1
           Rint_rel=PAle{klen}.Rint_rel;
          else
           Rint_rel=0;
          end
          if isfield(PAle{klen},'Rax')==1
           Rax=PAle{klen}.Rax;
          else
           Rax=1;
          end
          if Rax==1
           ShaM=7;
          else
           ShaM=4;
          end
          if ifp<=-11
%          'lens_sub prima in ass_parn', keyboard
%          'lens_sub prima in ass_parn', keyboard
%          'lens_sub prima in ass_parn', keyboard
          end


%          'lens_sub prima in ass_parn', keyboard
         if isfield(PAle{klen},'pos')==1
          gra_le.pos=PAle{klen}.pos;
          gra_le.thg=PAle{klen}.thg;
          gra_le.thga=PAle{klen}.thga;
          gra_le.d=PAle{klen}.d;
          gra_le.LA=PAle{klen}.LA;
          gra_le.lambda=lambda;
         else 
          gra_le=[];
         end
   
% determino dove si trova il reticolo efficace, se c'e
 
          if ifp==-10 
%           ' prima di Ass_lente primo', pausak
          end
  %' prima  KLEN', keyboard
          Ass_lente
         end  %klen
         
%  ' asesso distemo  KLEN', keyboard
          itrle=find(diff(ipudvf)>1);
%          keyboard
          if klen==1
           itrle=length(ipudv);
           pu_p=1:ipudv(1)-1;
           pu_p1=1:ipudv(1);
           pu_d=ipudv(itrle)+1:length(dv);
           dilu=Dilu{1};
           ailu=Ailu{1};
           bilu=Ailu{1}*Raxv(1);
           nilu=Nilu{1};
           pu_fst=Filu{1};
           ShaM=ShMv(1);
           lensas1
          else
           apu_p=1:ipudv(1)-1;
           apu_p1=1:ipudv(1);
           apu_d=ipudv(itrle)+1:ipudv(itrle+1)-1;
           dilu=Dilu{1};
           ailu=Ailu{1};

           bilu=Ailu{1}*Raxv(1);

           nilu=Nilu{1};
           pu_fst=Filu{1};
           bpu_d=ipudv(end):length(dv);
           bpu_d=ipudv(end)+1:length(dv);
           dilu2=Dilu{2};
           ailu2=Ailu{2};
           bilu2=Ailu{2}*Raxv(2);
           nilu2=Nilu{2};
           pu_fst2=Filu{2};
           ShaM=ShMv(2);
           lensas2
           shavet_le(ipudvf(1:itrle))=0;
          end
          
          
%  ' asesso distemo ', keyboard
%  ' asesso distemo ', keyboard

         if Gac~=0
          if ifp==-10 
           ' prima di Ass_lente, seconda volta ', pausak
           keyboard
          end
          fi_replace=find(real(nv0)<0);
          n_mesav=sqrt(gra_le.n_me);
          if length(fi_replace)>0
           nversa=nv0(fi_replace(1));
          else
           fi_replace=firet;
          end 
         end 
          
          Ass_lente
          if ifp==-10 
           ' dopo  Ass_lente, seconda volta ', pausak
           keyboard
          end          
          if isfield(gra_le,'n_me')==1
           nrefret=-sqrt(gra_le.n_me);
           if iretpiano==1
            nv0(fi_replace,1)=nrefret;
           else
            nv0(fi_replace)=nrefret;
           end
          end 
         

            if ifp==-10
             ' dopo  reass LENTI', 
%             keyboard
           end         
       end  % lente
       
    if Ps.ifpstop==1   
%     disp(' loo 2 fine '), keyboard
    end

   reassign

%  end

 end  %length(par_in)>0
 end
%' doppo lente', keyboard

if ~exist('gra_le') 
  gra_le.Gac=1;
end
if gra_le.Gac==0 & igraef_new==1
 ireturn=0;
 deter_gamma
end

if ifp==-10
%'dopo PARAMETRI', keyboard
end

if ifp==-10
'% ricalcola i campi in z'
%keyboard
end
global irec_fie iord_long num_long lp1c
if length(irec_fie)==0
 irec_fie=1;
end
if length(iord_long)==0
 iord_long=0;
end

 if  iord_long==0
%  fi_z
'entro fiz2', keyboard
  fiz_2
 else
%  fsr_sub
%'%  fsr_subnu  eliminato in Ass_parn'
 end

if ifp==-10
'dopo ricalcolo campo '
%keyboard
end

%disp(' fine asspar'), keyboard
%disp(' fine asspar'), keyboard

% if iset_MESA==1
%  mesa
% end
 
  if iset_MESA==1
    mesa_venCompleto
%    disp(' fine MESA'), keyboard
  end

% move 1
% Set parametri per autoconsistenza (rela_new.m)
%
% per struttura planare - autoconsistenza

  nstratid=nmir.b;
  Lvbr=Lpl.b.i;
  Luvb=Lpl.b.o;
  Lbb =Lpl.b.m;
  nvbr=npl.b.i;
  nuvb=npl.b.o;
  nbb =npl.b.m;
  rfd =npl.b.last;

  Lvb=[Lr.b; Lpl.b.i];
  nvb=[nr.b(:,1); npl.b.i];

  nstratiu=nmir.t;
  Lvtr=Lpl.t.i;
  Luv=Lpl.t.o;
  Lbt =Lpl.t.m;
  nvtr=npl.t.i;
  nuv=npl.t.o;

  nbt =npl.t.m;
  rfu =npl.t.last;

% per struttura non planare e/o  anisotropa

  Litn=Lr.t;
  nitn=nr.t.';
  global set_perd0
  if length(set_perd0)==0
   set_perd0=0;
  end
%  'iper', keyboard
  if set_perd0==1
   nip=nitn;
   fin=find(abs(imag(nip))<.1);
%   nip(fin)=nip(fin)-i*imag(nip(fin));
   nip(fin)=real(nip(fin));
   nitn=nip;
  end

%' qui perd', keyboard

  aitn=ar.t.';

  fmlsp=abs(frp.t);
  fmlst=abs(frp.t);
  fi1=find(fmlsp(:,1)==0);
  fmlst(fi1,1)=1;
  fi2=find(fmlsp(:,2)==1);
  fmlst(fi1,2)=0;



  Lib=Lr.b;
  nib=nr.b.';
  aib=ar.b.';
  if set_perd0==1
   nip=nib;
   fin=find(abs(imag(nip))<.1);
%   nip(fin)=nip(fin)-i*imag(nip(fin));
   nip(fin)=real(nip(fin));
   nib=nip;
  end
%  ' asspar', keyboard

  fmlsp=abs(frp.b);
  fmlsb=abs(frp.b);
  fi1=find(fmlsp(:,1)==0);
  fmlsb(fi1,1)=1;
  fi2=find(fmlsp(:,2)==1);
  fmlsb(fi1,2)=0;

  d=Lr.a*1e-6;
  L_QW=Lr.a;
  NQW=nmir.a;
  if NQW>1
   if length(Lr.t)>0
    dbar=Lr.t(length(Lr.t))*1e-6;
   else
    dbar=Lr.b(1)*1e-6;
   end
  else
   dbar=0;
  end
  niat=nr.a.';
  aiat=ar.a.';
  rqw=niat(1);
%'niat', keyboard  
  ipilat=0;
  if length(aiat)>1
   if aiat(2)>0
    ipilat=1;
   end
  end

% move 2

   pucavi=1;

   Litot=[Litn; d*1e6; Lib];
   aitot=[aitn.'; aiat'; aib.'];
   %nia
   nitot=[nitn.'; niat.'; nib.'];
   fmlstot=[fmlst; [1 0] ; fmlsb];
%   'asspar',   keyboard
%if ifp==-10 | Ps.ifpstop==1
 show_VCSELpro
%end

% per struttura completamente planare (perdite e acc. uscita in TL_ef.m)
%
% sopra
%
Lplit=Lpla.t.i;
nplit=npla.t.i;
%
Lplt=Lpla.t.m;
nplt=npla.t.m;
Nsplt=Nspla.t.m;
%
Lplut=Lpla.t.o;
nplut=npla.t.o;
%
%%
%
% sotto
%
Lplib=Lpla.b.i;
nplib=npla.b.i;
%
Lplb=Lpla.b.m;
nplb=npla.b.m;
Nsplb=Nspla.b.m;
%
Lplub=Lpla.b.o;
nplub=npla.b.o;

% save parpla
if ifp==-10
% 'salvato parpla', keyboard
end
if ifp==-11
%ifp=-10
end

 %'salvato parpla', keyboard
if  irec_fie==1 | lp1c==1
%if  irec_fie==1 
 if  iord_long==0
  fi_z
 else
%  fsr_sub_2010
%  fsr_sub_2010nu
%global lamm lp1 lp2 lp3 lastimm ilargeint Dlam0 demm iBEL Dlam_mod gastimm gazet gacen
global lamm lp1 lp2 lp3 lastimm ilargeint Dlam0 demm iBEL gastimm gazet gacen
global la1Ds ga1Ds 
global la1Dr ga1Dr 

  ifpBEL=ifp;

if Ps.ifpstop==1
 'campi 1D'%, keyboard
% ipcam=input(' vuoi plot campi? ')
 ipcam=1;

 if length(ipcam)==0
  ifp=-4;
 else
  if ipcam==0
     ifp=-4;
  end
 end 
else  
 ipcam=0;
end

%'cont iBEL', keyboard
% if length(iBEL)==0
  if isfield(Ps,'iBEL')==1
   iBEL=Ps.iBEL;
  else
   'manca iBEL in Ass_parn', keyboard
  end
% end 
%  iBEL=Ps.iBEL;
  if iBEL==-1
%   lambda=lamm{lp1,lp2}(1);
   fsr_sub_2010BEL
   lambda=.945;
%   lambda=.946;
   %gth1d=1000;
   lastimm(lp1,lp2)=lambda;
   
%   ' ver -1', keyboard
  elseif iBEL==1963

%   lambda=lambda0*1e6;  
%   lambda=1.5579;   
   %' ad hoc pri ' , keyboard   
   fsr_sub_2010BEL
   
   if ga1Ds<ga1Dr
    lambda=la1Ds;
   else
    lambda=la1Dr;
   end
   lastimm(lp1,lp2)=lambda;
   ' ad hoc' , keyboard   
%   ' ver -1', keyboard
  elseif iBEL==1995

%   lambda=lambda0*1e6;  
%   lambda=1.5579;   
%   ' nuovi campi CMM ' , keyboard   
	lambda=lambda0*1e6;

	iret_BW=0;
%ifp=-10

   fisiz=find(aitot>0);
   adu=sort(aitot(fisiz));
   fi=find(diff(adu)>0);
   if length(adu)>0 & length(fi)==0
    alat=adu(1);
   else
   alat=[];
   for kl=1:length(fi)
    alat=[alat; adu(fi(kl)+[0 1])];
   end 
   end
      fi=find(diff(alat)>0);
      alat=alat([fi; length(alat)]);
      
 ' nuovi campi CMM ' , keyboard        

	[iem,gain, lamr,zet,Ez,nz,NQW_ef,par_grat,Neff_ro,nzu,gatot]=campCMMI(lambda,dv,nv,xm,radii,fst,ifp,shavet,iauto,anyf,iret_BW,L_i,n_i,alat,ipolar,rr);
  ' nuovi campi CMM ' , keyboard   
 
 if T~=0
 
   fiT=find(zet>=zdis(1) & zet<=zdis(end));
   dz=[0 diff(zet)];
   Ezedz=abs(Ez).^2.*dz;   
   noE=sum(Ezedz);   
   Iz=sum(Ezedz(fiT));   
   confzT=Iz/noE;

   Nrot=T(:,1)';
   %dnT=Nrot*confzT;
   
   dnT=Nrot*confzT;
   
   roeff=ro_inT;
   dN_att=del_n_ag*gatot*ianti_gui;

%'passo lat', keyboard
   clear nlat
   fi0=find(roeff<aiat(1));
   fi=find(roeff<alat(1));
   nlat(fi)=Neff_ro(1);
   for kl=1:length(alat)-1
    fi=find(roeff>=alat(kl) & roeff<alat(kl+1));
    nlat(fi)=Neff_ro(kl+1);
%    figure, plot(roeff(fi),nlat(fi)), pausak
   end
   fi=find(roeff>=alat(end));
   nlat(fi)=Neff_ro(end);   
if ifp==-10 | Ps.ifpstop==1   
figure, 
nag=nlat;
nag(fi0)=nag(fi0)+dN_att;
subplot(211), plot(roeff,nlat,roeff,nag,roeff,sqrt(nlat.^2+dnT),'linewidth',2),
xlabel(' radial coord (um)')
ylabel(' Weighed index ')
k0cm=-2*pi/lambda*1e4;
subplot(212), plot(roeff,imag(nlat)*k0cm,roeff,imag(sqrt(nlat.^2+dnT))*k0cm,'linewidth',2), 
xlabel(' radial coord (um)')
ylabel(' Weighed Losses (1/cm) ')
pausak
figure, 
plot(roeff,nlat,'linewidth',2),
xlabel(' radial coord (um)')
ylabel(' Weighed index ')
pausak

%figure, plot(zet, nzu-3.5,zet,4*abs(Ez).^2,zet, nzu), keyboard
 figure,
 subplot(212)
 plot(zet,nzu(:,:),zet,3.5*abs(Ez).^2)
xlabel(' Long coord (um)')
 subplot(211)
 plot(zet,nzu(:,:),zet,3.5*abs(Ez).^2)
 axis([zdis(end)-.1 zedis(end)+.7 0 4])

title(' Standing Wave & Ref. index profiles')
pausak

figure, 
plot(zet,abs(Ez).^2,zet(fiT),abs(Ez(:,fiT)).^2,'r.')
title(' Temperature affecting ')
xlabel(' Long coord (um)')

	Ezc=Ez;
	zc=zet;
	nc=nz(:,1);

%  	save sa Ezc zc nc

keyboard
end

nlatT=nlat+dnT;
if ifp==-10
figure, plot(roeff,nlat,roeff,nlatT,'linewidth',2), 
   ' doponuovi campi CMM ' , keyboard  
end

global Ther   
Ther.roeff=roeff;
Ther.nlat=nlat;
Ther.dnT=dnT;

end

if ifp==-10 | Ps.ifpstop==1
 figure,
 plot(zet,nzu(:,:),zet,3.5*abs(Ez).^2)
 xlabel(' Long coord (um)')
 ylabel(' Standing Wave & Ref. index profiles')
 title([' Gain= ',num2str(gain,3),' (1/cm)         Wav= ',num2str(lamr*1000,3),' (nm)'])
 pausak
% 'qui per campo', keyboard
end
   
   %if ga1Ds<ga1Dr
   % lambda=la1Ds;
   %else
   % lambda=la1Dr;
   %end
   lambda=lamr;
   lastimm(lp1,lp2)=lamr;
   gastimm(lp1,lp2)=gain;
   gth1d=gain;
   uL(2)=NQW_ef;
   uL(1)=1;
   %' ad hoc' , keyboard   
   %' ver -1', keyboard

  elseif iBEL==2012

%   lambda=lambda0*1e6;  
%   lambda=1.5579;   

%   ' iBEL 2012 ' , keyboard  
unpack_0



  lam0=lambda0*1e6;
  kt=0;
   

      if n_i(end)~=rfd
       L_i=[L_i; 0];
       n_i=[n_i; rfd];
      end   
   
   L_imic=L_i;
   

   
   dlam=.03;   % fraction of lam0 where res is searched
%   ' 1D', keyboard
if ifp==-10
   ' 1D', keyboard
end
   rr1d=1;
   [gi,la,zet,nz,Ez,Hz,uL]=CMM_1D(lam0,L_imic,n_i,iat,icav,fiQW,ifp,dlam,rr1d,kt,ibast,par_grat,ipcam);

   lastimm(lp1,lp2)=la;
   gth1d=gi;
   lambda=la;
   la1D=lambda;
   Dlam_mod=Dlam0;
   NQW_ef=uL(2)/uL(1);
%   uL(2)=NQW;
%   uL(1)=1;
   %' ad hoc' , keyboard   
%   ' ver -1', keyboard
  elseif iBEL==0
   fsr_sub_2010nuu
   'qui iBEL = 0', keyboard  
  elseif iBEL==1
%  'qui', keyboard

   fsr_sub_2010BEL
   la1D=lambda;
   li=la1Ds;
   gi=ga1Ds;
   if length(li)==3
    lambda=li(2);   
   else
     lambda=li;   
   end
   lastimm(lp1,lp2)=lambda;
   gastimm(lp1,lp2)=gth1d;

   ilargeint(lp1,lp2)=iBEL;   
   Dlam_mod=Dlam0;   
   laintv=lambda*1000-Dlam_mod(1:2)   
   %'ver BEL', keyboard

  elseif iBEL==2
     fsr_sub_2010BEL
     Dlam_mod=Dlam0;   
     li=la1Ds;
     gi=ga1Ds;
     if lp1==1
      lambda=li(2);   
     else
      lambda=lamm{lp1-1,lp2}(1); 
      iv=[1 3];
      [du,im]=min(gi(iv));
      dl=-diff(li([2 iv(im)]))*1000; 
      Dlam_mod(1:2)=Dlam0(6)*[-1 1]*dl;
      Dlam_mod(3)=Dlam0(5);
     end
     lastimm(lp1,lp2)=lambda;
     gastimm(lp1,lp2)=gth1d;
  
     ilargeint(lp1,lp2)=iBEL;   

     laintv=lambda*1000-Dlam_mod(1:2) 
%    'dopo BEL 2 ', keyboard
   
%    'qui BEL 1', keyboard
  elseif iBEL==10000
%  'qui', keyboard

   Dlam_mod=Dlam0;
   fsr_sub_2010BEL
   la1D=lambda;
   ilargi=1
   lastimm(lp1,lp2)=lambda;
   ilargeint(lp1,lp2)=ilargi;
   
    'qui BEL 10', keyboard



  elseif iBEL==3
   fsr_sub_2010BEL
   lambda=.888
   Dlam_mod(1)=0
   ifp_post=1   
   'dopo BEL 3 ', keyboard
  elseif iBEL==4
   fsr_sub_2010BEL
   Dlam_mod(1)=-10
   Dlam_mod(2)=10
%   'dopo BEL 3 ', keyboard
  elseif iBEL==10 | iBEL==11
  
   fsr_sub_2010BEL
   li=la1Ds;
   gi=ga1Ds;
   iv=[1 3];
   [du,im]=min(gi(iv));
   dl=abs(diff(li([2 iv(im)]))); 
   lambda=li(2);   

  % dl=diff(li([2 iv(im)])); 
   
   
   PE1=Dlam0(1)*1000;
   PE2=Dlam0(2)*1000;
   L1=dl*PE1;
   L2=dl*PE2;   
   Dlam_mod=Dlam0;  
   Dlam_mod(1:2)=sort([L1 L2]);  
%'Dlam_mod', keyboard
   
   if iBEL==11
     if lp1>1   
       lambda=lamm{lp1-1,lp2}(1);
       Dlam_mod=Dlam0;
       Dlam_mod(3)=Dlam_mod(5);
       Dlam_mod(1:2)=Dlam_mod(6)*[-1 1]*dl*1000;
     end 
   end
   

   lastimm(lp1,lp2)=lambda;
   gastimm(lp1,lp2)=gth1d;

   ilargeint(lp1,lp2)=iBEL;   
   laintv=lambda*1000-Dlam_mod(1:2)
%' iBEL 10', keyboard

  elseif iBEL==24
  
   fsr_sub_2010BEL
   li=la1Ds;
   gi=ga1Ds;
   iv=[1 3];
   [du,im]=min(gi(iv));
   dl=-diff(li([2 iv(im)])); 
   lambda=li(2);   
   dl=-min(diff(li));
  % dl=diff(li([2 iv(im)])); 
   
   
   PE1=Dlam0(1)*1000;
   PE2=Dlam0(2)*1000;
   L1=dl*PE1;
   L2=dl*PE2;   
   Dlam_mod(1:2)=sort([L1 L2]);  
   

   lastimm(lp1,lp2)=lambda;
   gastimm(lp1,lp2)=gth1d;

   ilargeint(lp1,lp2)=-1;   
   laintv=lambda*1000-Dlam_mod(1:2)
 
     elseif iBEL==30
     
      fsr_sub_2010BEL
      li=la1Ds;
      gi=ga1Ds;
      iv=[1 3];
      [du,im]=min(gi(iv));
      dl=-diff(li([2 iv(im)])); 
      [gdu,isog]=sort(gi);   

      lambdavet=li(isog(1:2));   
      lambda=li(isog(1));   
      lambdac=li(2);   
      dl=-min(diff(li));
     % dl=diff(li([2 iv(im)])); 
      
      
      PE1=Dlam0(1)*1000;
      PE2=Dlam0(2)*1000;
      L1=dl*PE1;
      L2=dl*PE2;   
      Dlam_mod(1:2)=sort([L1 L2]);  
      
   
      lastimm(lp1,lp2)=lambdac;
      gastimm(lp1,lp2)=gth1d;
   
      ilargeint(lp1,lp2)=iBEL;   
   laintv=lambdac*1000-Dlam_mod(1:2)
%' iBEL 30', keyboard
     elseif iBEL==31 | iBEL==32 
     
      fsr_sub_2010BEL
      li=la1Ds;
      gi=ga1Ds;
      iv=[1 3];
      [du,im]=min(gi(iv));
      dl=-diff(li([2 iv(im)])); 
      [gdu,isog]=sort(gi);   

      lambdavet=li(isog(1:2));
      if iBEL==31
       lambdavet=li(isog(1));   
      end 
      lambda=li(isog(1));   
      lambdac=li(2);   
      dl=-min(diff(li));
%      dl=diff(li([2 iv(im)])); 
      
      
      PE1=Dlam0(1)*1000;
      PE2=Dlam0(2)*1000;
      L1=dl*PE1;
      L2=dl*PE2;   
      Dlam_mod(1:2)=sort([L1 L2]);  
      
     if lp1>1   
       lambdac=lamm{lp1-1,lp2}(1);
       Dlam_mod=Dlam0;
       Dlam_mod(3)=Dlam_mod(5);
       Dlam_mod(1:2)=Dlam_mod(6)*[-1 1];
       lambdavet=lambdavet(1);
     end 
      lastimm(lp1,lp2)=lambdac;
      gastimm(lp1,lp2)=gth1d;
   
      ilargeint(lp1,lp2)=iBEL;   
      
      laintv=lambdac*1000-Dlam_mod(1:2)
%' iBEL 31', keyboard


  elseif iBEL==20
  
   fsr_sub_2010BEL
   li=la1Ds;
   gi=ga1Ds;

%   li=la1Dr;
%   gi=ga1Dr;
   
   iv=[1 3];
   [du,im]=min(gi(iv));
   dl=-diff(li([2 iv(im)])); 
  % dl=diff(li([2 iv(im)])); 
   
   
   PE1=Dlam0(1)*1000;
   PE2=Dlam0(2)*1000;
   L1=dl*PE1;
   L2=dl*PE2;   
   Dlam_mod(1:2)=sort([L1 L2]);  
   
   lambda=li(2);
   lastimm(lp1,lp2)=lambda;
   gastimm(lp1,lp2)=gth1d;

   ilargeint(lp1,lp2)=-1;   
   laintv=lambda*1000-Dlam_mod(1:2)
%   'dopo BEL 20 ', keyboard
  elseif iBEL==21
   fsr_sub_2010BEL
   li=la1Ds;
   gi=ga1Ds;

%   li=la1Dr;
%   gi=ga1Dr;
   
   iv=[1 3];
   [du,im]=min(gi(iv));
   dl=-diff(li([2 iv(im)])); 
   
   PE1=Dlam0(1)*1000;
   PE2=Dlam0(2)*1000;
   L1=dl*PE1;
   L2=dl*PE2;   
   Dlam_mod(1:2)=sort([L1 L2]);  
   
   lambda=li(2);

   gastimm(lp1,lp2)=gth1d;

   ilargeint(lp1,lp2)=iBEL;   
   laintv=lambda*1000-Dlam_mod(1:2)
  
      if lp1>1
       lambda=lamm{lp1-1,lp2}(1);
       Dlam_mod=Dlam0;
       Dlam_mod(3)=Dlam_mod(5);
       Dlam_mod(1:2)=Dlam_mod(6)*[-1 1];
      end
   lastimm(lp1,lp2)=lambda;
%    'iBEL = 21', keyboard
  elseif iBEL==-10
   fsr_sub_2010BEL

   gastimm(lp1,lp2)=gth1d;
       Dlam_mod=Dlam0;

   lambda=Dlam_mod(7);
   lastimm(lp1,lp2)=lambda;
  elseif iBEL==-11
   fsr_sub_2010BEL

   gastimm(lp1,lp2)=gth1d;
   Dlam_mod=Dlam0;
   if lp1==1  
    lambda=Dlam_mod(7);
   else
    lambda=lamm{lp1-1,lp2}(1);
   end
    Dlam_mod(1:2)=Dlam_mod(6)*[-1 1];    
    Dlam_mod(3)=Dlam_mod(5);    
    lastimm(lp1,lp2)=lambda;   
          freq=0;
          TL_ef
          W=0;
    %      save Par1D uL z_Temp E_Temp Lef nmean W fTras fTrasinf fPinf fPsup Perd St_wave
          global uL z_Temp E_Temp Lef nmean W fTras fTrasinf fPinf fPsup Perd St_wave
      'TL', keyboard

  elseif iBEL==100
  
%   'iBEL=100', keyboard
   lambda=lambda0*1e6;
   fsr_sub_2010BEL
   lastimm(lp1,lp2)=lambda;
   gastimm(lp1,lp2)=gth1d;
   gazet(lp1,lp2)=uL(2);
   gacen(lp1,lp2)=uL(1);
%   'qie', keyboard
% ' load stop '
%   load stop
%   if istop==1
%    keyboard
%   end   
     
 elseif iBEL==101
    
 ' load stop '
  load stop
  if istop==1
   keyboard
  end      
     NP_sopra=Ps.NP_sopra;
     NP_sotto=Ps.NP_sotto;    
    
%     'iBEL=101', keyboard
     lambda=lambda0*1e6;
     L_pair=lambda/real(mean(n_i))*1000;
     
     par_grat.itetm=1;
     
     fsr_sub_2010BEL
     lastimm(lp1,lp2)=lambda;
     gastimm(lp1,lp2)=gth1d;
     la1Dr=lambda;
     ga1Dr=gth1d;
     
     if dret>0
      par_grat.itetm=2;
      fsr_sub_2010BEL
      lastimm1(lp1,lp2)=lambda;
      gastimm1(lp1,lp2)=gth1d;
     end 
      la1Ds=lambda;
      ga1Ds=gth1d;  
%      Litot(1:2)
%      [la1Dr la1Ds]
%      [ga1Dr ga1Ds]
      
%      pausak
    save Par1D uL
  
  elseif iBEL==102 | iBEL==99 | iBEL==103
   if ifp==-10
   'calcolo la polarizzazione per Veronique'
   ['iBEL= ',num2str(iBEL)], keyboard

   end
     NP_sopra=Ps.NP_sopra;
     NP_sotto=Ps.NP_sotto;       
   fish=find(shavet==-8);
   if iBEL==103
      fish=find(shavet==-4);
	  if length(fish)>0
	  shavet(fish)=4;
	  end	  
	  if iLP==1
	  	  shavet(fish)=1;
		  shav.t(fish-1)=1;
	  end
   end
   figra=fish;
   if length(figra)==0
    clear par_grat
   end
   if isfield(Ps,'iBWnew')==1   
    iBWnew=Ps.iBWnew;   
   else
     iBWnew=0;
   end
   
   fia=find(anyf~=0);

   
   nv3D=nv;
   if length(figra)==1
    fireBW=figra;
	riA=radii.array{figra};
	clear par_grat
	if length(radii.array{figra})==13
    par_grat=radii.array{figra}{13};
	else
	par_grat.period=riA{5};
	par_grat.DC=riA{6}/riA{5};
	end
    ngrating=nv(fish,:);
    g_n1=ngrating(1);
    g_n2=ngrating(2);
    par_BW=par_grat;
    Period=par_BW.period*1000;
    lambda1000=lambda*1000;
    if length(lambda1000)==0
     lambda1000=lambda0*1e9;
    end
    thickness=Litot(fireBW-1)*1000;
    DC=par_BW.DC;

%     'IM PARGREA ', keyboard
	 
    if fireBW<fiCav(1)
     nin=nitot(fireBW,1);
     nout=nitot(fireBW-2,1);
    else
     nout=nitot(fireBW,1);
     nin=nitot(fireBW-2,1);  
    end
    NModes=11;
     if ipolar==1
      npol=2;
     elseif ipolar==-1
      npol=1;
	 else
      npol=2; 	 
     end
     %'VETR', keyboard
    [neff,nBW,flagt] = f_EffectiveIndex(Period,lambda1000,DC,nin,nout,g_n1,g_n2,thickness,NModes);
    if iBWnew==0
	 Pf.ANYG=nBW;
     e_pa=nBW(1)^2;
     e_ve=nBW(2)^2;
     ng=nBW(npol);
    else
     e_pa=neff(1)^2;
     e_ve=neff(2)^2;  
	 Pf.ANYG=neff;
     ng=neff(npol);
    end    


     nv(fish,1)=ng;
	 FISH=fish-1;
     ng,     %keyboard
%     ng,     keyboard
%     ng,     keyboard
%     ng,     keyboard
 

   end   
   
%   'ver gara', keyboard
   
   if abs(ipolar)~=1
    polvet=[-1 1];
   else
    polvet=ipolar;
   end
   if iBEL==103
    polvet=polvet(1);
   end
 ipinc=0;   
 nvsave=nv;
 for ipolar=polvet
  ipinc=ipinc+1;   
 
%    par_grat.itetm=-ipolar;   

%    if ipolar==1 & length(fish)>0
   if ipinc==1
    if iLP==0 & iany==3 
%     if igraef_new==0
     'modifica 2018'
     anit=any.t;
     fiani=find(anit~=0);
	 iMa=imag(anit(fiani));
	 if length(find(iMa>1))>0
      nv(fiani+1,1)=imag(anit(fiani));
	 end 
%      nv(fiani+2,1)=anit(fiani);
%     end
    end

	if iBEL==103
	 anyf=anyf*0;
	 anyf(fish)=real(Pf.ANYG(1))+j*real(Pf.ANYG(2));
	 any.t=anyf(put);	
	 any.b=anyf(pub);	
	end
%    'ver anyf', keyboard    
   end 	
%    'ver lambda', keyboard    
   %if ~exist(
   lambda=lambda0*1e6;
   fsr_sub_2010BEL
   ipinchi=1;
   if ipolar==1
    ipinchi=2;
   end
   lambdavet(ipinc)=lambda2(ipinchi);
%    'ver gara', keyboard     
   gavet(ipinc)=G(ipinchi);
   uLv(ipinc,:)=uL;
  end 
   global gpla
   lastimm(lp1,lp2,lp3)=lambdavet(1);
   gastimm(lp1,lp2,lp3)=gavet(1);
   gpla=gth1d;
   gazet(lp1,lp2,lp3)=uLv(1,2);
   gacen(lp1,lp2,lp3)=uLv(1,1);
%    'qie', keyboard
%   nv=nvsave;

   freq=0;
   W=0;
   TL_ef   
   
   fatqw=uL(2)/uL(1)/NQW;
   if iBEL==99
    lambda=Ps.lambda_Temper/1000;
    lambdavet=ones(size(lambdavet))*Ps.lambda_Temper/1000;   
%     'verifico 99'
%    keyboard
   end 

%  ' load stop '
% 
%  ' load stop '
%   load stop
%   if istop==1
%    keyboard
%   end     
  
  end % iBEL

 end
 if ~exist('lambdavet')
  lambdavet=lambda;
 end
 if ~exist('dla_fre')
  dla_fre=0;
 end
% 'lambda =',
%'dopo Bel' , keyboard
if ~exist('i1D')
 i1D=0;
end 

if i1D==1
  eval(['save ' nomeloELM, ' lambda -append -v6' ]);
end  
end 

if exist('ifpBEL')
ifp=ifpBEL;
end
if i1D==1
 return
end
% disp(' fine ass_parn'), keyboard

shavet0=shavet;

%'fine Assparn', keyboard