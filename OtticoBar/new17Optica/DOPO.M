

gra=length(aou)-3;
co=polyfit(fou,aou,gra);
pfo=linspace(fou(1),fou(end),100);
af=polyval(co,pfo);
zed=roots(co);
df=diff(fou(1:2));
fia=find(imag(zed)==0 & (real(zed)<fou(end) & real(zed)>fou(1)));
if length(fia)==0
 fia=find(imag(zed)==0 & (real(zed)<fou(end)+df & real(zed)>fou(1)-df ));
end
zev=zed(fia);

% figure, plot(fou,aou,'g',ze,0,'wo',pfo,af,'r',zev,zev*0,'w*'), pausak

co=polyfit(fou,log10(gou),gra);
gf=10.^(polyval(co,pfo));
g0v=10.^(polyval(co,zev));

% figure, plot(fou,gou,'c',ze,gg0,'wo',pfo,gf,'m',zev,g0v,'w*'), pausak

for ks=1:length(zev)
 zel=zev(ks);
 [du,fi]=sort(abs(fou-zel));
 pu=[1:3];
 grl=length(pu)-1;
 fi=sort(fi(pu));
 fz=fou(fi);
 az=aou(fi);
 co=polyfit(fz,az,grl);
 zed=roots(co);
 fia=find(imag(zed)==0 & (real(zed)<fz(end) & real(zed)>fz(1)));
 zeu(ks)=zed(fia);
 co=polyfit(fz,log10(gou(fi)),grl);
 g0k=10.^(polyval(co,zed(fia)));
 g0u(ks)=g0k;
end
[g0s,ifs]=min(g0u);
zes=zeu(ifs);

ze=zes;
gg0=g0s;

% figure, plot(fou,aou,'r',zeu,g0u*0,'wo',zes,0,'w*'), pausak
% figure, plot(fou,gou,'c',zeu,g0u,'wo',zes,g0s,'w*'), pausak
