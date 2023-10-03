function [dilu,ailu,nilu,fstu,nver,gra]=lens_sub(th,h,r,N,nl,updown,ifp,Rel,Hre,Npair,Rflat,Rm_ring,gra)
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
%' in sub', keyboard
if length(Npair)==1
 Npair(2)=0;
end 

  Np_ad=Npair(2);
  Npair=Npair(1);

%  [dilu,ailu,nilu,fstu]=lens_subold(th,h,r,N,nl,updown,ifp,Rel,Hre,Npair,Rflat,Rm_ring,gra);
%  [dilu,ailu,nilu,fstu]=lens_vecchio(th,h,r,N,nl,updown,ifp,Rel,Hre,Npair,Rflat,Rm_ring)

%'% lens_suba fa tutti i casi, adesso  '
%'%IN LENSUB: prima di lens_suba fa tutti i casi, adesso  '
%keyboard
  %[dilu,ailu,nilu,fstu,nver,gra]=lens_suba(th,h,r,N,nl,updown,ifp,Rel,Hre,Npair,Rflat,Rm_ring,Np_ad,gra);
  [dilu,ailu,nilu,fstu,nver,gra]=lens_suba2016(th,h,r,N,nl,updown,ifp,Rel,Hre,Npair,Rflat,Rm_ring,Np_ad,gra);
% lens_suba fa tutti i casi, adesso

%' fine lensuba', keyboard
%' fine lensuba', keyboard

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

