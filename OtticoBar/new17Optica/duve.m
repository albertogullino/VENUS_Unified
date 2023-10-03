load
sgimp=1;
igint=0;
ifalso=-1;


% parte vecchia
%X=pdi;
%D_stri=aloc/kcav0;
%d_stri=bloc/kcav0;
%r_add=arig;
%
%if D_stri<0
% D_stri=-D_stri;
% orgra=pi/2;
%else
% orgra=0;
%end

%if d_stri<0
% d_stri=-d_stri;
% shiret=-D_stri/2;
%else
% shiret=0;
%end

%'grat'
%keyboard
% parte nuova
X=P.Rx;
Y=P.Ry;
Shap=P.shape;
D_stri=P.D;
d_stri=P.d;
r_add=P.Radd;
if P.orien~0
 orgra=pi/2;
 Xlo=X;
 Ylo=Y;
else
 orgra=0;
 Xlo=Y;
 Ylo=X;
end

if P.shif~0
 shiret=-D_stri/2;
else
 shiret=0;
end

   switch Shap
    case 'circle'
     Xma=Y;
    case {'square','rectangle'}
     Xma=sqrt(Y^2+X^2);
   end
npr=101*fix(Xma/D_stri);
rv=linspace(0,Xma,npr);

if ifp>-3 | ifp==-10 | ifp>=0

 ns=round(2*Xlo/D_stri);
 if shiret==0
  if is_even(ns)==1
   ns=ns+1;
  end
 else
  if is_even(ns)==0
   ns=ns+1;
  end
 end
  figure
%  yv0=[-d_stri -d_stri d_stri d_stri]/2;
  yv0=[-d_stri*ones(1,100) d_stri*ones(1,100)]/2;
  for k=1:ns
   Ym=(k-1-fix(ns/2))*D_stri-shiret;
   yvd=yv0+Ym;
   switch Shap
    case 'circle'
     dif=Xlo^2-(Ym)^2;
    case {'square','rectangle'}
     dif=Ylo^2;
   end
   if dif>0
    Xm=sqrt(dif);
    xpar=linspace(-Xm,Xm,100);
    xvd=[xpar fliplr(xpar)];
    if orgra==0
     yv=yvd;
     xv=xvd;
    else
     xv=yvd;
     yv=xvd;
    end
    if r_add>0
%    keyboard
     fiad=find(xv.^2+yv.^2<r_add^2);
     xv(fiad)=NaN;
     yv(fiad)=NaN;
    end
    plot3(xv,yv,ones(size(xv)),'w'), view(2), hold on
   end
  end
    if r_add>0
     fi=linspace(0,2*pi,length(xv));
     cirv=r_add*exp(j*fi);
     plot3(real(cirv),imag(cirv),ones(size(xv)),'w'), view(2)
    end

%  npr=1001;

  hold on,
   switch Shap
    case 'circle'
     polar(linspace(0,2*pi,100),max(rv)*ones(1,100),'r')
    case {'square','rectangle'}
     xplo=[X X -X -X X];
     yplo=[-Y Y Y -Y -Y];
     plot(xplo,yplo,'r')
   end
  axis equal
  drawnow
  if ifp>-3
   pausak
  end
end

muv=[0:2:2*nubesu];

    if r_add>0
     firv=find(rv>r_add);
     firv0=1:firv(1)-1;
    else
     firv=2:npr;
     firv0=1;
    end

AB=ones(length(rv),length(muv))*NaN;
im=0;
for mu=muv

im=im+1;

 if mu==0

  fat=pi/2;
  if r_add==0
   if shiret==0
    fat1=fat;
   else
    fat1=0;
   end
  else
    fat1=fat;
  end

  F_sh(firv0)=fat1;
  firv=100;
  for k=firv
   r=rv(k)
   a0=d_stri/(2*r);
   nstr=fix((r+d_stri/2-shiret)/D_stri);
%   nstr=fix((r+shiret)/D_stri);
   pausak
   som=0;
   som0=0;
   if  shiret==0
    if nstr==0 & shiret==0
     if a0<=1
      as0=asin(a0);
     else
      as0=fat;
     end
     som0=as0;
    else
     som0=asin(a0);
    end
   end
   for kj=1:nstr
    a1=(kj*D_stri+d_stri/2+shiret)/r;
    if a1<=1
     as1=asin(a1);
    else
     as1=fat;
    end
    a2=(kj*D_stri-d_stri/2+shiret)/r;
    as2=asin(a2);
    som=som+(as1-as2);
   end

   F_sh(k)=som+som0;
   if imag(F_sh(k))~=0
    pausak
   end
  end  %r
  if ifp>0
   ver=4*(F_sh*rv')*diff(rv(1:2))/(pi*max(rv)^2)/(d_stri/D_stri)
   pausak
  end

 else

  fat=0;
  F_sh(1)=fat;
  if r_add>0
   F_sh(firv0)=pi/4;
  end
  for k=firv
   r=rv(k);
   a0=d_stri/(2*r);
   nstr=fix((r+d_stri/2-shiret)/D_stri);
   som=0;
   som0=0;
   if shiret==0
    if nstr==0
     if a0<=1
      as0=asin(a0);
     else
      as0=fat;
     end
     som0=sin(mu*(as0+orgra))/mu;
    else
     as0=asin(a0);
     som0=sin(mu*(as0+orgra))/mu;
    end
   end

   for kj=1:nstr
    a1=(kj*D_stri+d_stri/2+shiret)/r;
    if a1<=1
     as1=asin(a1);
    else
     as1=fat;
    end
    a2=(kj*D_stri-d_stri/2+shiret)/r;
    as2=asin(a2);
    som=som+(sin(mu*(as1+orgra))-sin(mu*(as2+orgra)))/mu;
   end
   F_sh(k)=som+som0;
  end  %r
  F_sh=F_sh;
 end
 AB(:,im)=2*F_sh';

% figure, plot(rv,F,'.g')
% pausak
end

r=rv;
if ifp>1
 if length(muv)>6
  figure, plot(r,AB(:,1:6)), hold on, plot(r,AB(:,7:length(muv)),'.-'),
 else
  figure, plot(r,AB),
 end
 title(' y0, m2, c4, r6, g8, b10, y12, m14, c16, r18, g20, b22')
 pausak
end
r=rv*kcav0;
rv=r;

  if ifp>10
   figure
  end
 for imu=pimu
  jmu=imu-meun;
  mu=mbv(imu);
  for inu=pimu
   jnu=inu-meun;
   nu=mbv(inu);
   if (nu+mu)/2-fix((nu+mu)/2)==0
    dmn=abs(mu-nu);
    fim=find(dmn==muv);
    if length(fim)==1
     mfatd=AB(:,fim);
     A(:,jmu,jnu)=AB(:,fim);
    else
     disp('errore A mu in sha_grat ')
     pausak
    end
    dmn=abs(mu+nu);
    fim=find(dmn==muv);
    if length(fim)==1
     mfats=AB(:,fim);
     B(:,jmu,jnu)=AB(:,fim);
    else
     disp('errore B mu in sha_grat ')
     pausak
    end
   end  %if
  end
 end

%a=asav;
%disp(' fine sha_grat '), keyboard
