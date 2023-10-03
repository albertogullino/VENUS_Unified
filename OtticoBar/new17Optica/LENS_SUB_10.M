function [dilu,ailu,nilu,fstu,nver]=lens_sub(th,h,r,N,nl,updown,ifp,Rel,Hre,Npair,Rflat,Rm_ring,gra,PAlens)
%'in lensub', keyboard
if nargin<=12 
 gra=[];
end
if length(gra)==0
 gra.pos=0;
 gra.thg=0;
 gra.thga=0;
 gra.d=0;
 gra.LA=0;
end
%lens_suba fa tutti i casi, adesso
%' in sub', keyboard
%lens_suba fa tutti i casi, adesso
if length(Npair)==1
'%lens_suba fa tutti i casi, adesso'
 [dilu,ailu,nilu,fstu]=lens_subold(th,h,r,N,nl,updown,ifp,Rel,Hre,Npair,Rflat,Rm_ring,PAlens);
 %'esco lensubold', keyboard
 nver=0;
else 
  Npair(2)=0;
  Np_ad=Npair(2);
  Npair=Npair(1);
  [dilu,ailu,nilu,fstu,nver]=lens_suba(th,h,r,N,nl,updown,ifp,Rel,Hre,Npair,Rflat,Rm_ring,Np_ad,gra);
% lens_suba fa tutti i casi, adesso

%' fine lensuba', keyboard
end

fiva=find(fstu(:,2)>0);
fstu=fstu(fiva,:);
ailu=ailu(fiva,:);
dilu=dilu(fiva,:);
nilu=nilu(fiva,:);
per0=sum(th(1:2))/1000;

itro=find(diff(fstu(:,2))~=0);
itro=[itro; length(fstu)];
ki=1;
xd=[];
yd=[];
dsu=cumsum(dilu)/1000;
sup=0;
for k=1:length(itro)
 ku=itro(k);
 pu=ki:ku;
 ydu=dsu(pu);
 xdu=ailu(pu,:);
 RP=fstu(ki,2);
 if RP==0
  RP=1;
 end
 for kk=1:RP
  if kk>1
   sup=sup+per0;
  end
  xd=[xd; xdu];
  yd=[yd; ydu+sup];
 end
% if kk<RP
%  sup=sup+per01;
%  xd=[xd; xdu];
%  yd=[yd; ydu+sup];
% end
 ki=ku+1;
end
if ifp==-10
 figure
 plot(xd,yd,'.'), 
 grid, axis equal
 pausak
 %figure, plot(fstu(:,2),'.')
 %pausak

end

return
%vecchio modo
%' in sub', keyboard
 if Npair(2)==0 
   ' NON ha reticolo 1', keyboard
  [dilu,ailu,nilu,fstu]=lens_subold(th,h,r,N,nl,updown,ifp,Rel,Hre,Npair,Rflat,Rm_ring,gra);
 else
  Np_ad=Npair(2);
  Npair=Npair(1);
  [dilu,ailu,nilu,fstu]=lens_suba(th,h,r,N,nl,updown,ifp,Rel,Hre,Npair,Rflat,Rm_ring,Np_ad,gra);
 end 

