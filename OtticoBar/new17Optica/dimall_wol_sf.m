M_max=1.7;  %150
%M_max=1.5;  %200
M_max=1.5;  %250
M_max=2.5;  %fo
N_max=6
clear Gv Mv Lv Dv
idisper=0;
iplo=idisper
iff=0;
imod_M2=1
ifig=imod_M2;
le=length(demm);
Pp=PPlot{1};
XP=Pp.XP;
YP=Pp.YP;
X=Pp.X;
Y=Pp.Y;
xri=xro;
xdx2=xri.^3*diff(xro(1:2));
xdx=xri.*diff(xro(1:2));
X=PPlot{1}.X;
z=10e3;
tet=X(:,1)/180*pi;
xf=z*tan(tet)';
dx=diff(xf);
dx=[dx(1) dx];
xdxf=xf.*dx;
xdx2f=xf.^3.*dx;

iy=1;
isez=2;
fadox=3;
ica=isez;

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
fi_x=input(' punto ? (return tutti)  ')

for iy=1:size(demm,2)
 if length(fi_x)==0
  fi_x=find(isnan(Dp(:,iy))==0);
 end 
for ix=fi_x'
% 'qui dimo_subw', keyboard

dimo_subw
La1(ix,:)=La1Dr{ix};

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

%xdx2=xri.^3*diff(xro(1:2));
%xdx=xri.*diff(xro(1:2));
I=abs(Ex0.^2);
% si2=20^2
% I=exp(-2*xro'.^2/si2);
Es=xdx2*I;
En=xdx*I;
sig2=Es./En;
%'m2', keyboard

%z=10e3;
%tet=X(:,1)/180*pi;
%xf=z*tan(tet)';
%dx=diff(xf);
%dx=[dx(1) dx];
%xdx=xf.*dx;
%xdx2=xf.^3.*dx;
Esf=xdx2f*abs(Ef0.^2);
Enf=xdxf*abs(Ef0.^2);
sig2f=Esf./Enf;

dsi0=2*sqrt(sig2*2);
dsi0f=2*sqrt(sig2f*2);

M2=pi/(4*lambda*z)*dsi0.*dsi0f;
%'m2', keyboard
if ifig==1
f3=figure,
set(f3,'pos',[119    28   348   406])
f1=figure,
set(f1,'pos', [ 140   196   339   526])
end
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
laver=1000*lambda-dlav;
fiv=find(M2<M_max & M2>1);
%[du,fivd]=sort(gainv);
%%[du,fivd]=sort(M2.*gainv);
%fiv=fivd(1:N_max);
%'Me', keyboard

M2i=M2(fiv);
gainvi=gainv(fiv);
delfi=delf(fiv);
lavi=laver(fiv);
%'qui M2', keyboard

if iff==1
figure, plot(X(:,1),Ef0*exp(2),'linewidth',2), grid, title(' FF profiles / e^2 ')
xlabel(' Angle ')
axis([0 10 0 8])
pausak
end

cor='ymcgrbw';

%for k=1:le
kc=ix
if ifig==1
figure(f3)
subplot(211), plot(laver,gainv,'o',lavi,gainvi,'.')
title([' K =',num2str(kc)])
subplot(212), plot(laver,M2,'o',lavi,M2i,'.')
pausak

figure(f1)
subplot(311), plot(Dpp(kc),gainv,'o'), ylabel('Thres. Gain')
subplot(312), plot(Dpp(kc),dlav,'o'), ylabel('Detuning')
xlabel(YL)
%punt=1:length(Dpp);
subplot(313), plot(Dpp(kc),M2,'o'), ylabel('M2')
xlabel('# in vector ')
%' qui prima campo', keyboard
for kpl=puk
f2=figure;
set(f2,'pos', [479   255-30*(kpl-1)   675   462])
subplot(211),
 plot(abs(Aov{kpl}),['.-r']), hold on,
 plot(abs(Aqv{kpl}),['-g'])
 title([' Gain= ',num2str(gainv(kpl)),'     lam=  ',num2str(laver(kpl)),'     M2=  ',num2str(M2(kpl))])
 xlabel(' k index')
 ylabel(' Norm. Coefficients')
 subplot('position',[0.1 0.1 .3 .3])
% if ica==1


 titl='|EX|^2 Cavity';
 map_fnew(XP,YP,abs(Ecv(:,:,kpl)).^2,aax,Cug.x,Cug.y,Cug.z,titl,ibar,0)
 axis on
 axis([-1 1 -1 1]*dox*fadat)
 xlabel(' x (micron)')

% else
 subplot('position',[0.4 0.1 .3 .3])
 titl='|EX|^2 output';
 map_fnew(XP,YP,abs(Exv(:,:,kpl)).^2,aax,Cug.x,Cug.y,Cug.z,titl,ibar,0)
% end
 axis on
 axis([-1 1 -1 1]*dox*fadou)

 xlabel(' x (micron)')

 subplot('position',[0.65 0.13 .3 .3])
 titl='Far-Field';
  surf(X,Y,(Efv(:,:,kpl)).^2),
  shading('interp'), view(0,90),
  axis square, axis equal, grid,
  xlabel('angle (deg)')
 % pausak
 axis([-1 1 -1 1]*axli),
% pausak
end %ifig
 end
if ifig==1
 ' prima di chiudere', keyboard
 close(f1+1:f2)
end
 M2=M2';
 gainv=gainv';

% ' quo', keyboard

 [du,iso]=sort(delfi);
%iso=1:length(iso);
Ds=delfi(iso);
 Ms=M2i(iso);
 Gs=gainvi(iso);
 pu=1:length(iso);
 Mv(pu,ix)=Ms;
 Gv(pu,ix)=Gs;
 Dv(pu,ix)=Ds;
end
end
fi=find(Mv==0);
Mv(fi)=NaN;
Gv(fi)=NaN;
Dv(fi)=NaN;

fi=find(isnan(Dp)==0);
Dxx=Dp(fi);
Dx=Dp(fi)-Dp(1);
fi=find(Dxx==1000);
if length(fi)==1
 Dxx(fi)=-1;
 Dx=Dxx-Dxx(1);
end
%Dx=Dp-Dp(1);

h=figure, 
set(h,'pos',[120    20   361   732])
subplot(411), plot(Dx,Mv,'.-'), title(' M2 ')
subplot(412), plot(Dx,Gv,'.-'), title(' Gain')
wl=1000*lambda./(1+Dv);
subplot(413), plot(Dx,wl,'.-'), title(' Wavelenght')
subplot(414), plot(Dx,-diff(wl),'.-'), title(' WL spacing')
pausak

%return
Lvp=fliplr(1./(1+Dv')*1000*lambda);
Gvp=fliplr(Gv');
dloop=Dx;
'cont ', keyboard
alig_res3d