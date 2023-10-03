global igraef_new
%'deter gamma qui', keyboard
if igraef_new==0
   gra_le.Gac=0;
   gra_le.ru=0;
   Gac=0;
   return
end   
if ireturn==1
  return
end

%' quo', keyboard
if ~exist('iretpiano')
 iretpiano=0;
end 
 if  iretpiano==0
  firet=find(real(nv0(:,1))<0);
 else
  firet=find(shavet(:,1)==6);
 end
  if length(firet)==0
   firet=find(shavet(:,1)==6);
   iretpiano=1;
  end
  if length(firet)>0
   if firet(1)==2 | firet(end)==length(nv0)-1
    return
   end
  end 
  
  if length(firet)>0
   fiaut=find(iauto(:,1)==2);
   if firet(1)>fiaut
    iposret=1;  %sotto
   else
    iposret=-1;  %sopra
   end
% trovo gli indici per il gamma dopo il reticolo
   if iposret==-1
     fiind=2:firet(1)-1;
%     fiind=2:firet(1);
     fi_i=firet(end)+1;
     fi_u=1;
   else
     fiind=firet(end)+1:length(dv)-1;
     fi_i=firet(1)+1;
     fi_u=length(dv);     
   end

   if ireturn==1
      gra_le.ru=nv0(fi_u,1);
' prima di return', keyboard
' prima di return', keyboard
      
   end
   repet=fst(fiind,:);
   nvzo=nv0(fiind,:);
   dvzo=dv(fiind,:);
   idif=diff([0; repet(:,1)]);
   fie=find(idif~=0);
   if fie(1)>1
    isetto=[1; fie; length(repet)+1];
   else
     isetto=[fie; length(repet)+1];
   end
   pustra=[];
   kle=1;
   while kle<=length(repet)
    if repet(kle,1)==0
     pustra=[pustra;  kle]; 
     kle=kle+1;
    else
     nss0=repet(kle,1);
     nss=kle+[0:nss0-1]';
     nrep=repet(kle,2);
     kle=kle+nss0;
     pustra=[pustra;  repmat(nss,nrep,1)]; 
    end
%    fii=[isetto(kpu):isetto(kpu+1)-1]';
%    pustra=[pustra;  repmat(fii,repet(fii(1),2),1)]; 
   end   

%'ver pistar', keyboard   
   
%   for kpu=1:length(isetto)-1
%    fii=[isetto(kpu):isetto(kpu+1)-1]';
%    pustra=[pustra;  repmat(fii,repet(fii(1),2),1)]; 
%   end
   
   if iposret==-1
    pustra=flipud(pustra);
   end
   dc=dvzo(pustra)/1000;
   nc=nvzo(pustra,1);
   
   if ifp==-10
    ncdu=[nc nc];
    dcdu=[zeros(size(dc)) dc];
    ncd=reshape(ncdu.',prod(size(ncdu)),1);
    dcd=reshape(dcdu.',prod(size(ncdu)),1);
   
    figure, plot(cumsum(dcd),ncd), 
    title(' qui struttura di cui si calcola il coeff di riflessione per il calcolo del ret. equivalente')
    pausak
   end
   
   r_inc=nv0(fi_i,1);

%   r_inc=nv0(fi_i,1);
   r_outc=nv0(fi_u,1);

%   nci=r_inc;
%   dci=.1447;
%   [Gac]=gaemms(0,0,lambda,dci,nci,[],[],0,r_outc,rr,iLP,dc,nc,r_inc); 

%   [Gac]=gaemms(0,0,lambda,dc,nc,[],[],0,r_outc,rr,iLP,[],[],r_inc); 
%   [Gac]=gaemms(0,0,lambda,dc,nc,[],[],0,r_outc,rr,iLP,[],[],rr); 
if iretpiano==0
 r_incc=dunr(2);
else
 r_incc=nv0(firet,1);
end
if ifp==-10
%' locale per settem2', keyboard
end
   [Gac]=gaemms(0,0,lambda,dc,nc,[],[],0,r_outc,rr,iLP,[],[],r_incc); 
   gra_le.Gac=Gac;
   gra_le.ru=r_inc;
%   gra=gra_le;
   if ifp==-10
   ' cont QTR', keyboard
   end
  end
   if iretpiano==1
       ary=radii.array{firet};

        Pgra.D=ary{5};
        Pgra.d=ary{6};
        Pgra.Ry=ary{7};
        Pgra.Rx=ary{8};
        Pgra.shape=ary{9};
        Pgra.orien=ary{2};
        Pgra.shif=ary{3};
        Pgra.Radd=ary{4};
        if length(ary)>=10
         Pgra.grap=ary{10};
        else
         Pgra.grap=0;
        end
        if length(ary)>=11
         Pgra.Rext=ary{11};
        else
         Pgra.Rext=0;
        end
%        tt=Pgra.D;
%        gra.LA=Pgra.D;
%        gra.DC=Pgra.d/Pgra.D;
%        t1=Pgra.d;
%        ngra=nv0(firet,:);
        thick=dv(firet)/1000;
        n1=nv0(firet,1);
%        n2=nv0(firet+1,1);
        n2=nv0(firet,2);
        nu=n2;
        period=Pgra.D;
        t1=Pgra.d;
        t2=period-t1;
        DC=t1/period;

%       ' fine quio', keyboard

%        period=tt;
%        DC=t1/tt;
%        n1=ngra(1);
%        n2=ngra(2);
      
   
      ifpgi=ifp;
      ifpgi=-4;
      itetm=3;
      Nmodi=31;
%      'le_gra_inf giro_eq', keyboard 
      [n_palim,n_velim,n_pa,n_ve]=giro_eq(lambda,thick,period,DC,n1,n2,ifpgi,n1,nu,Gac,iretpiano,itetm,Nmodi);
     if ifp==-10
      'le_gra_planar', keyboard
     end
     n_me=((n_ve^2+n_pa^2)/2);
     n_di=((n_pa^2-n_ve^2)/2); 
     gra_le.n_me=n_me;
     gra_le.n_di=n_di;
     ex=n_ve^2;
     ey=n_pa^2;
  end