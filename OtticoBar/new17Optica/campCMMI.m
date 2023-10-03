function [iem,gain, lamr,zet,Ez,nz,NQW_ef,par_grat,Neff_ro,nzu,gaqw]=campCMMI(lambda,dv,nv,xm,radii,fst,ifp,shavet,iauto,anyf,iret_BW,L_i,n_i,alat,ipol,rr);

%ifp=-10
%rr=n_i(1);
NPlam=30;    % punti scansione lambda
dlam=1e-1;
dlam=5e-2;
kt=0;

%'in cmm', keyboard

[shu,Lu,nu,fiQ,fiCav,fiQW,atot,ntot,iate]=unpack_fun1(dv,nv,xm,radii,fst,ifp,shavet,iauto,anyf);

nv(fiQ,:)=real(nv(fiQ,:));
dpac=dv;
npac=nv(:,1);
fpac=fst;

%'contr Fie', keyboard
%unpack_1


  fish=find(shavet(:,1)==6);
  ibast=fish;
  ibaste=find(shu==6);
  i2pol=ipol;
  if length(fish)>0
   i2pol=1;
  ngrating=nv(fish,:);
  par_grat.r1=ngrating(1);
  par_grat.r2=ngrating(2);
  par=radii.array{fish};
  n1g=ngrating(1);
  n2g=ngrating(2);
  period=par{5};
  DC=par{6}/period;
  t1=period*DC;
  t2=period-t1;
  par_grat.per=period;
  par_grat.DC=DC;
  NModi=21;
  par_grat.NModi=NModi;
  par_grat.iret_BW=iret_BW;
  par_grat.r_in=(nv(ibast-1,1));
  par_grat.r_out=(nv(ibast+1,1));
 else
  par_grat=[];
 end

%L_inm=L_i*1000;
%n_inm=n_i;
%
%if length(Lu)-1==length(L_inm)
% L_inm=[L_inm; Lu(end)];
% n_inm=[n_inm; nu(end)];
%else
%'errore struttura in Campi_CMM', keyboard
%end

L_inm=Lu;
n_inm=nu;
n_inm=ntot;
a_inm=atot;

kt=0;  %modi direzione normale strati

%'cont', keyboard
%ifp=-10

if ipol==1 | ipol==2

  itetm=2;
  par_grat.itetm=itetm;
%[gam,lam,zetm,Ezm,nzm,gatotm,NQW_efm,azm]=fiezCMMn(lambda,L_inm,n_inm,fiQ,fiCav,fiQW,ifp,dlam,NPlam,rr,kt,ibast,par_grat,a_inm);
[gam,lam,zetm,Ezm,nzm,gatotm,NQW_efm,azm]=fiezCMMIn(lambda,L_inm,n_inm,iate,dpac,npac,fpac,fiQ,fiQW,ifp,dlam,NPlam,rr,kt,ibast,ibaste,par_grat,a_inm,fiCav);
%  [gam,lam]=fiezCMMn(lambda,L_i,n_i,fiQW,fiCav,fiQ,ifp,dlam,NPlam,rr,kt,ibast,par_grat,ato);

NQW_efm
%'NQW', pausak

end
if ipol==-1 | ipol==2

  itetm=1;
  par_grat.itetm=itetm;
 [gae,lae,zete,Eze,nze,gatote,NQW_efe,aze]=fiezCMMIn(lambda,L_inm,n_inm,iate,dpac,npac,fpac,fiQ,fiQW,ifp,dlam,NPlam,rr,kt,ibast,ibaste,par_grat,a_inm,fiCav);

end

if ipol==1
 gae=gam;
 lae=lam;
 zete=zetm;
 Eze=Ezm;
 nze=nzm;
 gatote=gatotm;
 NQW_efe=NQW_efm;
 aze=azm;
elseif ipol==-1
 gam=gae;
 lam=lae;
 zetm=zete;
 Ezm=Eze;
 nzm=nze;
 gatotm=gatote;
 NQW_efm=NQW_efe;
 azm=aze;
end

if abs(gae)<abs(gam)
 iem=1;
 gain=gae;
 lamr=lae;
 Ez=Eze;
 nz=nze;
 zet=zete;
 NQW_ef=NQW_efe;
 az=aze;
 gaqw=gatote;
else 
 iem=2;
 gain=gam;
 lamr=lam;
 Ez=Ezm;
 nz=nzm;
 zet=zetm;
 NQW_ef=NQW_efm;
 az=azm;
 gaqw=gatotm;
end 
 


if ifp==-10
figure, 
subplot(211)
plot(zete,abs(Eze).^2*3, zete, nze,'r'), 
 lab='TE';
 title(['CMM: ',lab,'  lambda_{res} = ',num2str(lae),' Gth TE= ',num2str(gae)]), 
 ylabel(' Intens/Index')

subplot(212)
 lab='TM';
 plot(zetm,abs(Ezm).^2*3, zetm, nzm,'r'), 
 title(['CMM: ',lab,'  lambda_{res} = ',num2str(lam),' Gth TM= ',num2str(gam)]),  
 xlabel(' Long. coord. (um)')
 ylabel(' Intens/Index')
  pausak
 end


%'fare nuovo nz', keyboard

nzu=[];
alati=[0 alat' 1000];
sl=length(alat)+1;
for k=1:length(az)
 azi=az(k,:);
 azid=azi;
 nzii=nz(k,:);
 nziis=nzii;
 if length(find(azi>0)>0)
  izer=0;
  ib=find(azi>0);
  ib0=find(azi==0);
  nzii(ib0)=0;
  fima=find(nzii~=0);
  nzii=nzii(fima);
  azi=[0 azi(ib) 1000];
 else
  izer=1;
 end 
 
 if izer==0
  nzi=[];
  fia=find(diff(azi)~=0);
  azi=azi(fia);
  nzii=nzii(fia);    
   for ks=1:length(azi)-1
    fi=find(alati>=azi(ks) & alati<azi(ks+1));
    if length(fi)>0
     nad=nzii(ks)*ones(size(fi));
     nzi=[nzi nad];
    end
   end
  if length(nzi)<length(alati)-1
   fiu=find(diff(azid)<0)+1;
   if isempty(fiu)==1
    fiu=length(nziis);
   end   
   nzi=[nzi nziis(fiu)*ones(1,(length(alati)-1-length(nzi)))];
  end 
 else 
   nzi=ones(1,sl)*nzii(1);
 end
 nzu(k,1:length(nzi))=nzi;
% 'ref tras', k
% pausak
end

dz=[0 diff(zete)];
Ezedz=abs(Eze).^2.*dz;
noE=sum(Ezedz);
Neff_ro1=sqrt(Ezedz*nz.^2/noE);
%for k=1:length(Neff_ro1)
%  rdu=nzu(:,k);
% d1=find(diff(rdu)~=0);
% fiva=d1(1):d1(end);
% Neff_rou(k)=(Ezedz(fiva)*rdu(fiva)/noE);
%end 
Neff_rou=(Ezedz*nzu/noE);

Neff_ro=Neff_rou;
%' fine Campi per eff. index ', keyboard
