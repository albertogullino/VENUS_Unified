clear isez
if ~exist('dox')
 dox=xro(end)*.8;
end
if dox==0
 dox=3;
end
path(path,'D:\USERS\PIERLUIGI\MVCSEL\new03')
ix=input(' ix  ');
if length(ix)==0
ix=1;
end
iy=input(' iy  ');
if length(iy)==0
iy=1;
end

if length(dox)>1
 doxa=dox(ix);
end 
isez=input(' sezione? [0 centro   return,output]  ');
if length(isez)==0
 isez=1;
end
%['        L air      D BHJ       D etch          R1         R2        Rc']
%[D1(ix,iy) D2(ix,iy) D7(ix,iy) D8(ix,iy) D9(ix,iy) D5(ix,iy) ], pausak
%[D8(ix,iy) D9(ix,iy) D10(ix,iy) D11(ix,iy)], pausak
ica=isez;
Pp=PPlot{ix,iy};
uL=Ppol.uL;
faGam=uL(2)/uL(1);
vg=3e10/rr*faGam;
vg=faGam;
%vg=3e10/3.3*7;

ldap=Pp.ldap;
ldap1=ldap+1;

iffz=0;


coc(1,1:3)='g.-';
coc(2,1:3)='r.-';
coc(3,1:3)='c.-';
coc(4,1:3)='m.-';
coc(5,1:3)='y.-';
coc(6,1:3)='b.-';
coc(7,1:3)='w.-';
sc=size(coc);
lco=sc(1);
%gainv=gamm{ix,iy};
dat=D2(ix,iy);
gaind=Pp.parmod;
delf=demm{ix,iy};
lg=length(delf);
if lg==1
 puk=1;
else
 puk=1:lg;
end
for k=puk
 gainv(k)=gaind{k}(1)/vg;
 dlav(k)=delf(k)*1000*lambda;
end
[du,ilas]=min(gainv);
pograp0=[50 50    258   546];

for ip=puk
Ao=Pp.Ap{ip};
Aq=Pp.ApQ{ip};
if ica==1
 Ex=Pp.E2xo{ip};
else
 Ex=Pp.E2xp{ip};
end
Ef=Pp.Ef{ip};
if iLP==0
 if ica==1
  Ey=Pp.E2yo{ip};
 else
  Ey=Pp.E2yp{ip};
 end
end
%XP=Pp.XP{ip};
%YP=Pp.YP{ip};
%X=Pp.X{ip};
%Y=Pp.Y{ip};
%Cug=Pp.Cug{ip};
XP=Pp.XP;
YP=Pp.YP;
X=Pp.X;
Y=Pp.Y;
Cug=Pp.Cug;
aax=D3(1,1)*1.5;
if isnan(aax)==1
 aax=5;
end
if aax==0
 aax=doxa;
end
ibar=1;
iaoff=0;
titl='';
hc=figure;
%if ip==1
%pograp=[50 50 300 450];
%else
%pograp=[350 50 300 450];
%end
pograp=pograp0;
pograp(1)=pograp(1)+(ip-1)*900/puk(end);
set(gcf,'Position',pograp);

figure(hc)
subplot(411)
if iffz==0
titl='FF';
 surf(X,Y,(Ef).^2),
 shading('interp'), view(0,90),
 axis square, axis equal, grid,
% pausak
 axli=40;
 axis([-1 1 -1 1]*axli/2),
 title(titl)
else
 if iLP==0
  titl='Ez';
  map_fnew(XP,YP,Ez,aax,Cug.x,Cug.y,Cug.z,titl,ibar,iaoff)
 end
end

%map_fnew(X,Y,Ef,[],[],[],[],titl,ibar,iaoff)
subplot(412)
titl=[' alim ', num2str(alim),' Nk= ',num2str(length(Aq)/2)];
map_fnew(XP,YP,Ex,aax,Cug.x,Cug.y,Cug.z,titl,ibar,iaoff)
axis([-1 1 -1 1]*doxa*1.5)
if iLP==0
subplot(413)
titl='EY';
map_fnew(XP,YP,Ey,aax,Cug.x,Cug.y,Cug.z,titl,ibar,iaoff)
axis([-1 1 -1 1]*doxa*1.5)
end
%axis([-10 10 -10 10])
%gain=gamm{ix,iy}(ip);
gain=gainv(ip);
dlam=dlav(ip);

%stri=[' wav = ',num2str(lambda/(1+ze),'%0.6e')];
%text(2,9,stri);
%if ixy==1
% strip='X-pol';
%elseif ixy==-1
% strip='Y-pol';
%elseif ixy==0
% strip='LP mode';
%end
%text(5,7,strip);
hold off
subplot(414)
plot(abs(Ao/max(Ao)),coc(ip-fix((ip-1)/lco)*lco,:)), hold on,
plot(abs(Aq/max(Aq)),'w-')
 puau=1:length(Ao);
 puini=ldap1(1:end-1);
 if puini(end)<puau(end)
  puaud=puau(puini);
  plot(puaud,abs(Ao(puini)/max(Ao)),'co'),
 end

axis([0 length(Ao) 0 1])
stri=[' gain = ',num2str(gain,'%0.4e'),' wav = ',num2str(dlam,'%0.4e')];
%text(-9,9,stri);
ht=text(-.2,1.1,stri);
if ip==ilas
 set(ht,'color','y')
end
pausak
end
%


% figure
% surf(X,Y,10*log10((Ef/Ef(1,1)).^2)),
% titl='FF in dB';
% colorbar
% shading('interp'), view(0,90),
% axis square, axis equal, grid,

% figure
% surf(X,Y,((Ef/Ef(1,1)).^2)),
% titl='FF ';
% colorbar
% shading('interp'), view(0,90),
% axis square, axis equal, grid,
