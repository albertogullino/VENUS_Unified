   [GGe2,GGm2,TTe2,TTm2]=gaperd(KK,freq,lambda,Lvbr,nvbr,...
                   Lbb,nbb,nstratid,rfd,rr,iLP,Luvb,nuvb,rr,fapes);
 Ge2=[];
 Gm2=[];
 Te2=[];
 Tm2=[];
 for imu=1:pasnu:nubes+1
   Ge2=[Ge2; GGe2];
   Gm2=[Gm2; GGm2];
   Te2=[Te2; TTe2];
   Tm2=[Tm2; TTm2];
 end

 Gad=[Ge2; Gm2];
 Trd=[Te2; Tm2];
                   