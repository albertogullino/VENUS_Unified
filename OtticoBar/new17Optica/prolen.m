%load lenu
%load lenn
%keyboard
if length(Npair)==1
% [dilu,ailu,nilu,fstu]=lens_subold(th,h,r,N,nl,updown,ifp,Rel,Hre,Npair,Rflat,Rm_ring);
% [dilu,ailu,nilu,fstu]=lens_vecchio(th,h,r,N,nl,updown,ifp,Rel,Hre,Npair,Rflat,Rm_ring)
 [dilu,ailu,nilu,fstu]=lens_add(th,h,r,N,nl,ifp,Npair);
else
 Np_ad=Npair(2);
 Npair=Npair(1);
 [dilu,ailu,nilu,fstu]=lens_suba(th,h,r,N,nl,updown,ifp,Rel,Hre,Npair,Rflat,Rm_ring,Np_ad);
end 
