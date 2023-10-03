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
faGam=Ppol.uL(2)/Ppol.uL(1);
vg=3e10/rr*faGam;
%vg=3e10/rr;

for iy=1:size(demm,2)
fi=find(isnan(Dp(:,iy))==0);
for ix=fi'
dimo_subw
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

xdx2=xri.^3*diff(xro(1:2));
xdx=xri.*diff(xro(1:2));
I=abs(Ex0.^2);
% si2=20^2
% I=exp(-2*xro'.^2/si2);
Es=xdx2*I;
En=xdx*I;
sig2=Es./En;
%'m2', keyboard

z=10e3;
tet=X(:,1)/180*pi;
xf=z*tan(tet)';
dx=diff(xf);
dx=[dx(1) dx];
xdx=xf.*dx;
xdx2=xf.^3.*dx;
Esf=xdx2*abs(Ef0.^2);
Enf=xdx*abs(Ef0.^2);
sig2f=Esf./Enf;

dsi0=2*sqrt(sig2*2);
dsi0f=2*sqrt(sig2f*2);

M2=pi/(4*lambda*z)*dsi0.*dsi0f;
%'m2', keyboard

f1=figure,
set(f1,'pos', [ 140   196   339   526])
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
subplot(311), plot(Dpp,gax,'linewidth',2), title(' Threshold ')
subplot(312), plot(Dpp,laver,'linewidth',2), title(' wavelength ')
subplot(313), plot(Dpp,aper_vet,'linewidth',2), title(' FF aperture deg ')
xlabel(' pillar variation (nm)')
pausak

figure, plot(X(:,1),Efcut,'linewidth',2), grid, title(' FF profiles / e^2 ')
xlabel(' Angle ')
axis([0 10 0 8])
pausak

cor='ymcgrbw';
f2=figure
set(f2,'pos', [479   255   675   462])
for k=1:le
figure(f1)
subplot(311), plot(Dpp,gax,'.-g',Dpp(k),gax(k),'ro'), ylabel('Thres. Gain')
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


 titl='|EX|^2 Cavity';
 map_fnew(XP,YP,abs(Ecv(:,:,k)).^2,aax,Cug.x,Cug.y,Cug.z,titl,ibar,0)
 axis on
 axis([-1 1 -1 1]*dox*fadat)
 xlabel(' x (micron)')

% else
 subplot('position',[0.4 0.1 .3 .3])
 titl='|EX|^2 output';
 map_fnew(XP,YP,abs(Exv(:,:,k)).^2,aax,Cug.x,Cug.y,Cug.z,titl,ibar,0)
% end
 axis on
 axis([-1 1 -1 1]*dox*fadou)

 xlabel(' x (micron)')

 subplot('position',[0.65 0.13 .3 .3])
 titl='Far-Field';
  surf(X,Y,(Efv(:,:,k)).^2),
  shading('interp'), view(0,90),
  axis square, axis equal, grid,
  xlabel('angle (deg)')
 % pausak
 axis([-1 1 -1 1]*axli/2),
 pausak
 clg
end

end

clear