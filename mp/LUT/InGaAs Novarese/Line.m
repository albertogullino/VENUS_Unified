
lambda=lambdavet(Ila)*1e-7;
Gc=G(:,Ila)';
Ec=Es(:,Ila)';

LaMat=(2*pi./lambda);    
gg=Gc./LaMat;  
nn=squeeze(Dep);
nn=nn(:,Ila)';
alfv=diff(nn)./diff(gg);
alfv=[alfv alfv(end)];

Gth=[990 1800];

Dv=1e-12*Densityv;
for kk=1:2
 [du,im]=min(abs(Gc-Gth(kk)));
 imf=im+[-2:2];
 imv(kk)=im;
 cog=polyfit(Dv(imf),Gc(imf),2);
 coe=polyfit(Dv(imf),Ec(imf),2);
 coa=polyfit(Dv(imf),alfv(imf),2);
 cog(end)=cog(end)-Gth(kk);
 nthv=roots(cog);
 [du,it]=min(abs(nthv-Dv(im)));
 nth(kk)=nthv(it);
 rsp(kk)=polyval(coe,nth(kk));
 alv(kk)=polyval(coa,nth(kk));
end

figure, 
subplot(121)
plot(Dv,Gc), hold on
O=ones(size(Dv))';
GV=O*Gth;
plot(Dv,GV,'k--')
plot(nth,zeros(size(nth)),'ro');
 rn=repmat(nth,2,1);
 L=[0 0; 2000 2000];
plot(rn,L,'r--')
grid
axis([0 3 0 2000])
xlabel(' N  (1e12/cm^2)')
ylabel('Gain')



vg=.8*1e10;
c=3e10;
Gam=.024;
eta0=[.5 .71];

hc=6.62e-34*c;
el=1.6e-19;  





Lin=1e-3*(Gam.^2*vg*hc)*Gth.*eta0./(4*pi*lambda).*rsp


subplot(122)
plot(Dv,Ec,nth,rsp,'ro')
xlabel(' N  (1e12/cm^2)')
ylabel('R_{sp}')
text(Dv(fix(imv(1)*1.2)),Ec(fix(imv(1))),['LPP ',num2str(Lin(2),2),' mWxMHz'])
text(Dv(fix(imv(2)*1.2)),Ec(imv(2)),['LPP ',num2str(Lin(1),1),' mWxMHz'])

text(Dv(fix(imv(1)*1.9)),Ec(fix(imv(1)*1.8)),['\alpha =  ',num2str(alv(2),2)])
text(Dv(fix(imv(1)*1.9)),Ec(fix(imv(1)*1.5)),['\alpha =  ',num2str(alv(1),1)])