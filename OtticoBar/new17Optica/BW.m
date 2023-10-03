  if igraef_new==0
%' sono qui in mat_mix', keyboard 
%' sono qui in mat_mix', keyboard 
   er1=ngra(1)^2;
   er2=ngra(2)^2;
   n_ve=tt*er1*er2/(t2*er1+t1*er2);
   n_pa=(t1*er1+t2*er2)/tt;
   n_me=((n_ve+n_pa)/2);
   n_di=((n_pa-n_ve)/2);
  else 
   n_me=gra_le.n_me;
   n_di=gra_le.n_di;
  end
  
  if P.orien==0
   Delta_eps=n_di; 
  else
   Delta_eps=-n_di; 
  end

%%                rintz=nitr(ianu)+seg_ret*ani_gr/2;  % eps_z = eps_x