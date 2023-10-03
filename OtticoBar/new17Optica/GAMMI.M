' gammi ',

if iFFsect==3
 if igamveb==1  %coeff. di riflessione sotto
   [GGe3,GGm3]=gaemms(KK,ze,lambda,Lvbr,nvbr,...
                   Lbb,nbb,nstratid,rfd,rr,iLP,Luvb,nuvb);
 else
  GGe3=ones(size(KK))*GGbext;
  GGm3=GGe3;
 end
else

 if igamveu==1  %coeff. di riflessione sopra
   [GGe3,GGm3]=gaemms(KK,ze,lambda,Lvtr,nvtr,...
                     Lbt,nbt,nstratiu,rfu,rr,iLP,Luv,nuv);

 else
  GGe3=ones(size(KK))*GGuext;
  GGm3=GGe3;
 end
end

 Ge3=[];
 Gm3=[];

 for imu=1:pasnu:nubes+1
   Ge3=[Ge3; GGe3];
   Gm3=[Gm3; GGm3];
 end

if iLP==0
 GGei=[Ge3; Gm3];
else
  if iem_like==0
   GGei=(Ge3+Gm3)/2;
  elseif iem_like==1
% TE
   GGei=Ge3;
% TM
  elseif iem_like==2
   GGei=Gm3;
  elseif iem_like==-1
   GGei=Ge3(1)*ones(size(Ge3));
  end
end
Gainti=GGei(Pusas);
