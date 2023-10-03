function [gazk,enan,fak,Tue,Tb,Gue,Gum,Gbe,Gbm,Lf,Lcav,ztot,Ez,Hz,indz,nmean,Perd,Ge,Gm,lQW,Ezm,KKie,KKim]=...
         th_scattu(fiQWi,fiCavi,L_i,n_i,rr,nb,nu,lambda0,fr,kv,iLP,ifp,iff,ibast,par_grat,NPXi)


global iscattu iorta 

%' iscattu', keyboard
if ~exist('NPXi')
NPXi=1;
end
enan=0;
if iscattu==0
[gazk,Iku,Ikb,fak,Tue,Tb,Gue,Gum,Gbe,Gbm,Lf,Lcav,ztot,Ez,Hz,indz,nmean,Perd,Ge,Gm,lQW,Ezm,KKie,KKim]=...
         th_scattu0(fiQWi,fiCavi,L_i,n_i,rr,nb,nu,lambda0,fr,kv,iLP,ifp,iff,ibast,par_grat,NPXi);
else         
 if iorta==1
% lambda0
%  'thscattu3'
  [gazk,enan,fak,Tue,Tb,Gue,Gum,Gbe,Gbm,Lf,Lcav,ztot,Ez,Hz,indz,nmean,Perd,Ge,Gm,lQW,Ezm,KKie,KKim]=...
         th_scattu3last(fiQWi,fiCavi,L_i,n_i,rr,nb,nu,lambda0,fr,kv,iLP,ifp,iff,ibast,par_grat,NPXi);    
 else
   'thscattu2'
 [gazk,enan,fak,Tue,Tb,Gue,Gum,Gbe,Gbm,Lf,Lcav,ztot,Ez,Hz,indz,nmean,Perd,Ge,Gm,lQW,Ezm,KKie,KKim]=...
         th_scattu2(fiQWi,fiCavi,L_i,n_i,rr,nb,nu,lambda0,fr,kv,iLP,ifp,iff,ibast,par_grat,NPXi);    
 end % iorta
end         
  