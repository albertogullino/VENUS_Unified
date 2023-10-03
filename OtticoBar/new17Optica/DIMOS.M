%path(path,'D:\USERS\PIERLUIGI\MVCSEL\new03')
%ix=input(' ix  ');
%iy=input(' iy  ');
icao=0;
     ix=kl;
     iy=ke;
%     ' s
Pp=PPlot{ix,iy};
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
gaind=Pp.parmod;
lg=length(gaind);
if lg==1
 puk=1;
else
 puk=1:lg;
end
for k=puk
 gainv(k)=gaind{k}(1);
end
[du,ilas]=min(gainv);
pograp0=[50 -50 300 750];

for ip=puk
Ao=Pp.Ap{ip};
Aq=Pp.ApQ{ip};
if icao==1
 Ex=Pp.E2xo{ip};
else
 Ex=Pp.E2xp{ip};
end
Ef=Pp.Ef{ip};
if iLP==0
 if icao==1
  Ey=Pp.E2yo{ip};
 else
  Ey=Pp.E2yp{ip};
 end
end

mod_chap
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
aax=4;
ibar=1;
iaoff=0;
titl='';
if idisp==1
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
titl='FF';
 surf(X,Y,(Ef).^2),
 shading('interp'), view(0,90),
 axis square, axis equal, grid,
 axli=60;
 axis([-1 1 -1 1]*axli/2),
 title(titl)

%map_fnew(X,Y,Ef,[],[],[],[],titl,ibar,iaoff)
subplot(412)
titl='EX';
map_fnew(XP,YP,Ex,aax,Cug.x,Cug.y,Cug.z,titl,ibar,iaoff)
if iLP==0
subplot(413)
titl='EY';
map_fnew(XP,YP,Ey,aax,Cug.x,Cug.y,Cug.z,titl,ibar,iaoff)
end
%axis([-10 10 -10 10])
%gain=gamm{ix,iy}(ip);
gain=gainv(ip);

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
axis([0 length(Ao) 0 1])
stri=[' gain = ',num2str(gain,'%0.4e')];
%text(-9,9,stri);
ht=text(0,0.6,stri);
if ip==ilas
 set(ht,'color','y')
end
end %idisp
end
