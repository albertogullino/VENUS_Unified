
Pp=PPlot{ix,iy};


if idisper==1
 P=Pp;
% save
% 'qui disper_sub', keyboard
kl=ix;
ke=iy;
% disper_sub
 disper_sub1
% clear
% load
 %pausak
end
ldap=Pp.ldap;
ldap1=ldap+1;
dat=D2(ix,iy);
gaind=Pp.parmod;
delf=demm{ix,iy};
gam=gamm{ix,iy};
lal=lamm{ix,iy};
lg=length(gam);
puk=[1:lg];
clear gainv dlav laverdu
for k=puk
 gainv(k)=gam(k)/vg;
 dlav(k)=delf(k)*1000*lambda;
 laverdu(k)=lal(k)*1000;
end
[du,ilas]=min(gainv);
pograp0=[50 -50 350 750];
[du,im]=min(gainv);
    gax(ix,iy)=gainv(im);
    del(ix,iy)=dlav(im);
    laver(ix,iy)=laverdu(im);
%    laverdu,    'laverdu', pausak

clear Ef0 Ex0 Exu Infv
for ip=puk
Ao=Pp.Ap{ip};
Aq=Pp.ApQ{ip};
%Aov(:,ix)=Ao/max(Ao);
%Aqv(:,ip)=Aq/max(Aq);
Aov{ip}=Ao/max(Ao);
Aqv{ip}=Aq/max(Aq);
if icx==1
 Ex=Pp.E2xo{ip};
 Ec=Pp.E2xp{ip};
else 
 Ex=Pp.E2yo{ip};
 Ec=Pp.E2yp{ip};
end 
 Ex=Pp.E2yp{ip};
 Ec=Pp.E2xp{ip};

Ef=Pp.Ef{ip};
if iLP==0
 if ica==1
%  Ey=Pp.E2yo{ip};
  Ey=Pp.E2xp{ip};
 else
  Ey=Pp.E2yp{ip};
 end
lfia=ceil(length(fian)/2);
lme=[1 lfia];
Efv(:,:,ip)=Ef;
Ef0du=mean(Ef(:,lme)');
Ef0(:,ip)=Ef0du/max(Ef0du);
Enfdu=mean((abs(Ex(:,lme)).^2+abs(Ey(:,lme)))');
Infv(:,ip)=Enfdu;
Ex0(:,ip)=Ex(:,1);
Ex0(:,ip)=Ec(:,1);
Exu(:,ip)=Ex(:,1);
Exv(:,:,ip)=Ex;
Ecv(:,:,ip)=Ec;
Eyv(:,:,ip)=Ey;
%'ins u', pausak
end

end

