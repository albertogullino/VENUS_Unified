le=length(demm);
iy=1;
isez=2;
fadox=3;
ica=isez;
iffz=0;
if ica==1
 fadox=3;
else 
 fadox=1.5;
end 
axli=5;
 fadat=1.3;
 fadou=3;

clear Efv Ef0 Ex0 Exv Ecv Eyv
coc(1,1:3)='g.-';
coc(2,1:3)='r.-';
coc(3,1:3)='c.-';
coc(4,1:3)='m.-';
coc(5,1:3)='y.-';
coc(6,1:3)='b.-';
coc(7,1:3)='w.-';
sc=size(coc);
lco=sc(1);
vg=3e10/rr;

for iy=1:size(demm,2)

for ix=1:le
dimo_sub
end

XP=Pp.XP;
YP=Pp.YP;
X=Pp.X;
Y=Pp.Y;
Cug=Pp.Cug;
aax=D3(1,1)*1.5;
if isnan(aax)==1
 aax=5;
end
ibar=1;
iaoff=0;

xri=xro;
%fi=find(xro>aax);
%xri(fi)=0;
xdx=xri.^2*diff(xro(1:2));
dx=ones(size(xri))*diff(xro(1:2));
Es=xdx*abs(Ex0.^2);
En=dx*abs(Ex0.^2);
sig2=Es./En;

z=10e3;
tet=X(:,1)/180*pi;
xf=z*tan(tet)';
dx=diff(xf);
dx=[dx(1) dx];
xdx=xf.^2.*dx;
Esf=xdx*abs(Ef0.^2);
Enf=dx*abs(Ef0.^2);
sig2f=Esf./Enf;

dsi0=4*sqrt(sig2);
dsi0f=4*sqrt(sig2f);

M2=pi/(4*lambda*z)*dsi0.*dsi0f;

f1=figure,
set(f1,'pos', [ 0    55   432   866])
%figure
%plot(abs(Aov),'.-'), hold on,
%plot(abs(Aqv),'-')
YLv{1}='spacer variation (nm)';
 YLv{5}='radius of curvature (micron)';
  YLv{13}=' Sphere hight (micron)';
 YLv{18}=' Sphere off-axis (micron)';
 YLv{20}=' out Pairs'; 
 
 Dpp=Dp;
if np1==1
 Dpp=Dp-Dp(1);
end
if np1==1 | np1==5
 YL=YLv{np1};
else
 YL='';
end
Dpp=Dpp(1:length(gax));
ff=figure;
set(ff,'pos',[612    57   672   881])
subplot(311), semilogy(Dpp,gax)
subplot(312), plot(Dpp,laver)
subplot(313), plot(Dpp,dsi0/2), title(' W0 ')
pausak

cor='ymcgrbw';
f2=figure
set(f2,'pos', [440   377   831   524])
for k=1:le
figure(f1)
subplot(311), semilogy(Dpp,gax,'.-g',Dpp(k),gax(k),'ro'), ylabel('Thres. Gain')
subplot(312), plot(Dpp,del,'.-y',Dpp(k),del(k),'ro'), ylabel('Detuning')
xlabel(YL)
punt=1:length(Dpp);
subplot(313), plot(punt,M2,'.-c',punt(k),M2(k),'ro'), ylabel('M2')
xlabel('# in vector ')
figure(f2)
subplot(211),
 plot(abs(Aov{k}),['.-r']), hold on,
 plot(abs(Aqv{k}),['-g'])
 title([' parameter= ',num2str(Dpp(k)),'          Point # ',num2str(k)])
 xlabel(' k index')
 ylabel(' Norm. Coefficients')
 subplot('position',[0.1 0.1 .3 .3])
% if ica==1


 titl='|EX|^2 Out';
 map_fnew(XP,YP,abs(Exv(:,:,k)).^2,aax,Cug.x,Cug.y,Cug.z,titl,ibar,0)
 axis on
 axis([-1 1 -1 1]*dox*fadat)
 xlabel(' x (micron)')

% else
 subplot('position',[0.4 0.1 .3 .3])
 titl='|EY|^2 Out';
 map_fnew(XP,YP,abs(Eyv(:,:,k)).^2,aax,Cug.x,Cug.y,Cug.z,titl,ibar,0)
% end
 axis on
 axis([-1 1 -1 1]*dox*fadat)

 xlabel(' x (micron)')

 subplot('position',[0.65 0.13 .3 .3])
 titl='Far-Field';
  surf(X,Y,(Efv(:,:,k)).^2),
  shading('interp'), view(0,90),
  axis square, axis equal, grid,
  xlabel('angle (deg)')
 % pausak
 axis([-1 1 -1 1]*axli*2),
 pausak
 clg
end

end

clear