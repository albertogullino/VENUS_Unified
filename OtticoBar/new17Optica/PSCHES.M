fi=find(abs(xu)<1e-5);
xu(fi)=0;
h0=h;
R=(r^2+h^2)/(2*h);
X=(R-h);
X0=X;
te=atan(r/X);
tex=pi/2-linspace(0,te,100)';
ci=R*exp(j*tex);
xc=real(ci);
yc0=imag(ci)-R;
for k=1:length(yi)
 yc(:,k)=yc0+yi(k);
end


figure,
plot(xu,yu,'.'), grid
pausak
%'cont', keyboard
%'cont', keyboard
%'cont', keyboard

sx=size(xu);
sx2=sx(2)-1;

xup=xu;
%fi=find(xu==r);
%xup(fi)=0;
ly=1;
le=1;
ks=0;
cle=0;

'ui', keyboard
while le~=0
xip=10000;
dux=[];
duy=[];
k=ks;
kf=0;
cle=cle+1;
 while k<sx(1)
 k=k+1;
% [cle k]
 fi=find(xup(k,:)>0);
 %pausak
 if length(fi)>0
  kf=kf+1;
 end
 if kf==1
  xua=xu(k,fi);
  yua=yu(k)*ones(1,length(fi));
  xup(k,fi)=0;
  ks=k;
 elseif kf>1
  if length(fi)>0
   xua=xu(k,fi(1));
   yua=yu(k);
  end
 else
 end
 %[min(xua) xip], pausak
  if kf>0
%   xua
%   yua
%   pausak
   if min(xua)<xip
    dux=[dux xua];
    duy=[duy yua];
    xip=min(dux);
    if kf>1
     xup(k,fi(1))=0;
    end
    le=length(find(xup~=0))
   else

    [xso,iso]=sort(dux);
    cux(1:length(dux),ly)=xso';
    cuy(1:length(dux),ly)=duy(iso)';
    ly=ly+1;
%    pausak
%    keyboard

    k=10000;
   end
  end

 end  %while k
end

%    [xso,iso]=sort(dux);
%    cux(1:length(dux),ly)=xso';
%    cuy(1:length(dux),ly)=duy(iso)';

cuys=cuy;
%cuy=cuys-cuy(end,1);
figure, plot(cux,cuy,'.')
pausak
s=size(cux);
figure
for k=1:s(2)
du=cux(:,k);
dy=cuy(:,k);
fi=find(du>0);
dux=du(fi);
duy=dy(fi);
dux1=[0; dux];
dux2=[dux; r];
du=[dux1 dux2];
dxt=reshape(du',prod(size(du)),1);
duy1=[duy(1); duy];
duy2=[duy; duy(end)];
duy1=[duy; duy(end)];
duy2=[duy; duy(end)];
dy=[duy1 duy2];
dyt=reshape(dy',prod(size(du)),1)+yu(1);
plot(dxt,dyt),
hold on
%'cot1', pausak
end
plot(xc,yc), pausak