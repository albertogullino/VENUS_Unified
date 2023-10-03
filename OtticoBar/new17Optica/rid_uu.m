global pMu0u
PREC_RID=0.05;
%PREC_RID=0.1;
%PREC_RID=0.02;
%PREC_RID=0.005;
PREC_RID=0.001;
if ipolar~=0
 PREC_RID=1e-10;
% PREC_RID=0.05;
else
 PREC_RID=20e-2;
% PREC_RID=10e-10;
end 
% PREC_RID=0;
%' entr rid_uu ', keyboard
%save sa
%keyboard
mbvero=nubesi:pasnu:nubesu;
mbdu=(nubesi:pasnu:nubesu)+.1;
if iLP==0 
 mbdu=[mbdu -mbdu];
end
if ipolar==0 
 mbdu=[mbdu mbdu*i];
end
if ifp~=-4, '% riduco la base', end
%keyboard
firid=find(shailoop==5);
if length(firid>0)
% irid_bas=1;
 if iLP==1
  fiLP=1;
  npam=1;
 else
  fiLP=2;
  npam=2;
 end
else
% irid_bas=0;
  fiLP=1;
  npam=1;
end

global ispeed
%' irid_bas ', keyboard
%irid_bas=0;
if irid_bas==0
 pre=.0;
else
 pre=PREC_RID;
end


 if length(irid_bas)==0
  irid_bas=1;
 end

% if ispeed==0 | ipolar==0
%  irid_bas=0;
% end
%'qui;' , keyboard
%if irid_bas==1
if irid_bas>=0
 %asce=axtot+aytot;
 asce=aytot;
 sa=size(asce);
 if min(sa)>1
%  [du1,ips1]=max(asce);
  [du1,ips1]=max(shavet);
  [du2,ips2]=max(du1);
  irim=ips2;
  icom=du1(ips2);
  if max(du1)==6
   [du1,ips1]=max(asce);
   [du2,ips2]=max(du1);
   irim=ips1(ips2);
   icom=ips2;   
  end
  %'vedo', keyboard
 else
  [du1,icom]=max(asce);
  irim=1;
 end
%' prima', keyboard 
if ~exist('cce')
 cce=0;
end 

if icom<4
if abs(cce)>0 & length(find(shailoop==4))>0
 icom=4;
 irim=1;
end
end 
 if prod(size(MKosp))>1
  Md1=MKosp{icom};
  Md2=MKosm{icom};
  if length(Md1)==2
   Md1=Md1{2};
   Md2=Md2{2};
  else 
   if iscell(Md1)
    Md1=Md1{1};
    Md2=Md2{1};   
   end 
  end
  sM=length(size(Md1));
  if sM==2
   irim=1;
  end
  if icom==11
   Kplot1=sum(abs(Md1(:,:,:)),3);
%   Kplot2=sum(abs(Md2(:,:,3),3));
  else 
   Kplot1=Md1(:,:,irim);
   Kplot2=Md2(:,:,irim);
  end 

  clear Md1 Md2
%  Kplot1=MKosp{icom};
%  Kplot2=MKosm{icom};
 else
  Md1=MKosp{1};
  Md2=MKosm{1};
  Kplot1=Md1(:,:,icom);
  Kplot2=Md2(:,:,icom);
  clear Md1 Md2
 end
 if ipolar==0
  
 end
 %' dopo', keyboard
 siK=size(Kplot1);
 if length(siK)>2
  Kplot1=reshape(Kplot1(:,:,irim),siK(1:2));
 end
 siK=size(Kplot2);
 if length(siK)>2
  Kplot2=reshape(Kplot2(:,:,irim),siK(1:2));
 end
iall=2; 
iallo=2; 
if ipolar~=0
 iall=1; 
end
if iLP==1
 iallo=1;
end
if pasnu==1
 if ipolar==1
  Dm=diag(Kplot1);
 elseif ipolar==-1
  Dm=diag(Kplot2);
 elseif ipolar==2
  Dm=diag(Kplot2)+diag(Kplot1);  
 elseif ipolar==0
% 'qui', keyboard
%  Dm=diag(Kplot1)+.1;
  pu=1:length(KK);
  if mm==0
   pu=pu+length(KK);
  end
  Dm=sum(abs(Kplot1))';
  fi=find(Dm<PREC_RID);
  Dm(fi)=0;
 end
 
 
else
 Dm=diag(Kplot1+Kplot2);
end
%  ' quit irim', keyboard
 if ifp>=1
  figure, plot(Dm),
 end

 lk=length(KK);
 pu=1:lk;
 Pus=[];
 ldap=[0];

 Dmv=Dm;
if iLP==10 & ipolar~=0
  Dm=(Dm(1:end/2)+Dm([1:end/2]+end/2))/2;
 end
 if irid_bas==0
  Dm=Dm+1;
 end
 Moe=[];
 for nm=1:npam:numodi*iallo*iall
  pui=pu+(nm-1)*lk;
  Du=abs(Dm(pui));
  if max(Du)>0
  Du=Du/max(Du);
  der=diff([0; Du]);
  fiu=find(Du-pre>0 & der>=0);
  pua=fiu(1):lk;
  if npam==1
   Pus=[Pus pui(pua)];
   ldap=[ldap ldap(end)+length(pua)];
   Moe=[Moe ones(1,length(pua))*mbdu(nm)];
  else
   Pus=[Pus pui(pua) pui(pua)+lk];
   ldad=ldap(end)+length(pua);
   ldap=[ldap ldad ldad+length(pua)];
  end
  end %if
 end
% ' fine ldap', keyboard
% if iLP==0
%  Pus0=Pus;
%  Pus=[Pus Pus+lk*numodi];
%  ldap=[ldap ldap+lk*numodi];
% end

 lknumodi=lk*numodi;
 lknumodi1=ldap(end);

 if iLP==0
  Pus0=Pus;
%  Pusas0=Pus;
 if iallo~=2
  Pus=[Pus Pus+lknumodi];
  ldapu=[ldap ldap(2:end)+lknumodi1];
  ldap=[ldap ldap+lknumodi1];
 else
  ldapu=ldap;
 end
 else
  ldapu=ldap;
 end

%'Pus settato ', keyboard
 if iredmat==1
  Kploto1=Kplot1;
  if ipolar==0
   Kploto2=Kplot1;
  else
   Kploto2=Kplot2;  
  end
  Kplot1=Kplot1(Pus,Pus);
  if ipolar==0
  Kplot2=Kplot1;
  else
  Kplot2=Kplot2(Pus,Pus);
  end
  Pusas=Pus;
  Pusa=Pus;
%  Pusasav=Pus;
  Pus=1:length(Pus);
  Pusasav=Pus;
  Pusasav=Pusas;
 else
  si2=si;
  Pusa=1:si2(1);
%  Pusa=1:si(1);
  Kploto1=Kplot1;
  Kploto2=Kplot2;
 end

Pusasf=Pus;

% ifp=1
 if ifp>=1
  lD=length(Dmv);
  putdu=1:lD;
  figure,
  plot(putdu,Dmv,putdu(Pus),Dmv(Pus),'r.')
  [length(Dmv) length(Pus)]
%  map(Kploto)
  map(Kplot1(Pus,Pus))
 end

else

 lKA=si2(1);
 Pus=1:lKA;
 numoa=0:numodi-1;
% ld=[numoa numodi+numoa]*nk1max;
% ldap=[ld ld+lKA];
 ld=[numoa numodi+numoa]*nk1max;
 ldap=[ld ld(end)+nk1max];

end

lPus=length(Pus);
if iredmat==0
 si21=si2(1);
else
 si21=lPus;
end

%Pusa=1:si21;
Pust=[Pus Pus+si21];

% set global IdeOo pMu0u pMu0u1 lKA nk1max nures sim0

lKAn=length(Pus);
 global ive
 ive=0;
 if iredmat==1
  ive=0;
 end
 if iany==0 | inuo_bas==0
  ive=1;
 end

%' iredmat', keyboard
if iredmat==0
 if iLP==0
  lKA=nk1max*numodi*2;
  nures=numodi*4;
  sim0=nk1max*2;
 else
  lKA=nk1max*numodi;
  nures=numodi*2;
  sim0=nk1max*2;
 end
  numoa=0:numodi-1;
  ld=[numoa numodi+numoa]*nk1max;
  ldu=[ld ld+lKA];
  pMu0d=[];
  for kdu=1:nures
   pMu0d=[pMu0d ldu+(kdu-1)*2*lKA*nk1max];
  end
 if ive==1
  pMu0u=pMu0d;
 else
  for ki=1:nk1max
   pMu0u{ki}=pMu0d+ki+(ki-1)*2*lKA;
  end
 end

  pMu0u1=[0 2*lKA^2; lKA lKA+2*lKA^2];

  ld1=1:nk1max;
  if iLP==0
   ldud0=[ld1];
   ldud=[ldud0 ldud0+lKA];
  else
   ldud0=[ld1 ld1+lKA];
   ldud=ldud0;
  end
  pMc=[];
  for kdu=ldud
   pMc=[pMc ldud+(kdu-1)*lKA*2];
  end
  ldus=ldud;

  ldud=ld1;
  pMei=[];
  for kdu=ldud
   pMei=[pMei ldud+(kdu-1)*lKA];
  end

else  %riduco matrici

 if iLP==0
  lKA=lKAn;
  nures=numodi*4;
  sim0=nk1max*2;
 else
  lKA=lKAn;
  nures=numodi*2;
  sim0=nk1max*2;
 end
  nk1maxv=diff(ldapu);
%  if iLP==0
   nk1maxv=[nk1maxv diff(ldapu)];
%  end
  numoa=0:numodi-1;
  ld=ldapu;
  ldu=[ld ld(2:end-1)+lKA];

%  pMu0u=[];
%  for kdu=1:nures
%   pMu0u=[pMu0u ldu+(kdu-1)*2*lKA*nk1maxv(kdu)];
%  end
  if iLP==0
   nktot=4*nk1max*numodi;
  else
   nktot=2*nk1max*numodi;
  end
  VPU=1:nktot^2;
  MPU=reshape(VPU,nktot,nktot);
  if ipolar~0
   Pusta=[Pusa Pusa+Pusa(end)];
  else 
   Pusta=Pusa;
  end 
  MPUR=MPU(Pusta,Pusta);

  numoa=0:numodi-1;
  ld=[numoa numodi+numoa]*nk1max;
  lKA1=nktot/2;

  ldu=[ld ld+lKA1];
    pMu0d=[];
    for kdu=1:nures
     pMu0d=[pMu0d ldu+(kdu-1)*2*lKA1*nk1max];
    end

  if ive==1
   pMu0u=pMu0d;
  else

ivemepu=0;
    for ki=1:nk1max
     pMu0du=pMu0d+ki+(ki-1)*2*lKA1;
%     tic
     MPUdu=MPU;
     MPUdu(pMu0du)=-MPU(pMu0du);
     MPURdu=MPUdu(Pusta,Pusta);
     fib=find(MPURdu<0);
%     toc
     if ivemepu==1
     tic
     Vdu=MPU(pMu0du);
     lV=length(Vdu);
     ic=0;
     clear lVd
     for kj=1:lV
      fi=find(Vdu(kj)==MPUR);
      if length(fi)==1
       ic=ic+1;
       lVd(ic)=fi;
      end
     end
     toc
     length(find(lVd-fib'~=0))
     [ki, length(lVd)], pausak
     else
      lVd=fib;
     end
     pMu0u{ki}=lVd;
    end
   end

  pMu0u1=[0 2*lKA^2; lKA lKA+2*lKA^2];

%  for
    imo=1;
    ld1p=ldapu(imo+[0 1])+[1 0];
    ld1=ld1p(1):ld1p(2);
    if iLP==0
     ldud0=[ld1];
     ldud=[ldud0 ldud0+lKA];
    else
     ldud0=[ld1 ld1+lKA];
     ldud=ldud0;
    end

    pMc=[];
    for kdu=ldud
     pMc=[pMc ldud+(kdu-1)*lKA*2];
    end

    ldus=ldud;
    ldud=ld1;
    pMei=[];
    for kdu=ldud
     pMei=[pMei ldud+(kdu-1)*lKA];
    end

end

 IdeOo=diag(ones(2*lKA,1));
 IdeOon=diag(ones(2*lKAn,1));
 si2=[si21,si21];

%if ifp~=-4

%end
%keyboard
%disp(' IdeOo')
%keyboard



%(JMO
%'wuo ap', keyboard
if ipolar==0 | iLP1==0
 KKvd=[KKv; KKv];
 ZEv2=[ZEv; ZEv];
 if iLP==0 
  KKvd=repmat(KKvd,iall,1);
  ZEv2=repmat(ZEv2,iall,1);
 end
else
 ZEv2=[ZEv];
 KKvd=KKv;
end
if ifp~=-4,
    if igr_app==0
        if iLP==1
            figure, plot(KKvd(Pusas),'.')
        else
            figure, plot(KKvd(Pus0),'.')
        end
    else
        figure, plot(KKvd,'.')
    end
end

if igr_app==0
    if iLP==1
        modeset=KKvd(Pusas);
    else
        modeset=KKvd(Pus0);
    end
else
    modeset=KKv;
end

disp(' RID_UU'),
[length(KKvd) length(Pus)]

%' fine rid_uu', keyboard
if ifp==-10
' fine rid_uu', keyboard
end
%pausak
%JMO)
