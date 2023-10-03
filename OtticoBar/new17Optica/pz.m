%load saz
load sazu
icum=0;
nzu=[];
alati=[0 alat' 1000];
sl=length(alat)+1;
for k=1:length(az)
 azid=az(k,:);
 azi=azid;
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
  %pausak
 else
  izer=1;
 end 
 
 if izer==0
%azi, pausak

  nzi=[];
 fia=find(diff(azi)~=0);
 azi=azi(fia);
 nzii=nzii(fia);  
   for ks=1:length(azi)-1
    fi=find(alati>=azi(ks) & alati<azi(ks+1))
    if length(fi)>0
      nad=nzii(ks)*ones(size(fi));
      nzi=[nzi nad];
    end
  if k==2775
  'cont', keyboard
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
  if k==2775
  'cont 1', keyboard
 end 
 nzu(k,1:length(nzi))=nzi;
 if icum==1
 %'ref tras', keyboard
 end
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