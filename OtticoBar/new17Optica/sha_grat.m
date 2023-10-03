%' inizio sha_grat', keyboard

if ~exist('igraef_new')
 if isfield(Ps,'igraef_new')==1
  igraef_new=Ps.igraef_new;
 else
  igraef_new=-1;
 end 
end

sgimp=1;
igint=0;
ifalso=-1;

%' in sha_grat', keyboard

global icambio_ref iNIexch

Pstrut=P;

clear A B AB F_sh

X=P.Rx;
Y=P.Ry;
Shap=P.shape;
D_stri=P.D;
d_stri=P.d;
r_add=P.Radd;
Rext=P.Rext;
if isempty(Rext)==1
 Rext=0;
end
if isstruct(Rext)==1
 Rext=0;
end
if Rext<0
 Rext=0;
end
if Rext>0
 iNIexch=1;
end

if P.orien~0
 orgra=pi/2;
 Xlo=X;
 Ylo=Y;
else
 orgra=0;
 Xlo=Y;
 Ylo=X;
end
i_grap=P.grap;
if igraef_new>=1
 i_grap=1;
 sha=0;
end 

global igr_app
if length(igr_app)==1
% if igr_app==1
%  i_grap=igr_app;
% end
end


if iNIexch==1
 nitot(ind_vet-1,1:2)=nitot(ind_vet-1,[2 1]);
 nv(ind_vet,1:2)=nv(ind_vet,[2 1]);
 nitn=nv(put,:).';
 nib =nv(pub,:).';
 d_stri=P.D-P.d;
end

%' grat ', keyboard
ireturn=0;
 if ~exist('ifirstgrat')
  ifirstgrat=0;
 end
if i_grap>=1
 if iany==0 & igr_app==1
  ifirstgrat=1;
 end
   %' SHAGRAT in', keyboard
 if ifp~=-4
%  ' SHAGRAT in', keyboard
 end
 if ~exist('t1')
  t1=P.d;
  t2=P.D-P.d;
  tt=P.D;
 end 
  rad=radii.a;
 er1=nv(ind_vet,1)^2;
 er2=nv(ind_vet,2)^2;
'% reticolo con denti paralleli a x'
%keyboard
if P.orien==1
 ey=tt*er1*er2/(t2*er1+t1*er2);
 ex=(t1*er1+t2*er2)/tt;
else
 ex=tt*er1*er2/(t2*er1+t1*er2);
 ey=(t1*er1+t2*er2)/tt;
end

 n_me=((ex+ey)/2);
 n_di=((ex-ey)/2);

%  ' SHAGRAT in', keyboard



%' qui grat', keyboard
if igraef_new>=1
%nuovo
  period=tt;
  DC=t1/tt;
  fi6=find(shavet==6);
  ngra=nv(fi6,:);
  n1=ngra(1);
  n2=ngra(2);
  thick=dv(fi6)/1000;

 if ifp==-10
  ' % determino dove si trova il reticolo efficace, in SHA_GRAT', pausak
  ' % determino dove si trova il reticolo efficace, in SHA_GRAT', pausak
  ' % determino dove si trova il reticolo efficace, in SHA_GRAT', pausak
 end
 iretpiano=1;
% ' sjagrat', keyboard
 ' sjagrat cambiata logica soluzione piana', 
 if igraef_new==2
  ' in sha_gra ritorno se Orta'
  return
 end
 if igraef_new==3
   nv(ind_vet,1)=sqrt(n_me);
   aysa=ay;
   fi2d  
  ay=aysa;
  nr1=ngra(1);
  nr2=ngra(2);
  %'gra if', keyboard
  thick=dv(fi6)/1000;
  return
 else
  nr1=ngra(1);
  nr2=ngra(2);
  thick=dv(fi6)/1000; 
   %'gra else', keyboard
%  return
 end
 
 
 if igraef_new==1
  deter_gamma  
  
%'dopo gamme', keyboard
 
%   if ifp==-10
%     'dopo deter_gamma in sha_grat: vedi Gac  '
%     pausak
%   end 
%   nu=gra_le.ru;
%   Gac=gra_le.Gac; 
%   if ifp==-10
%    ' prma di giroea fpremp', keyboard
%    ' prma di giroea fpremp', keyboard
%    ' prma di giroea fpremp', keyboard
%    ' prma di giroea fpremp', keyboard
%   end 
%
%   save Gac Gac n1 n2 thick lambda period
%
%   [n_palim,n_velim,n_pa1,n_ve1]=giro_eq(lambda,thick,period,DC,n1,n2,ifp,nu,Gac,iretpiano)  
%    ' prma di giroea fpremp dopo, A COSA SERVE ?', keyboard   
%    ' prma di giroea fpremp dopo, A COSA SERVE', keyboard   
%    ' prma di giroea fpremp dopo, A COSA SERVE', keyboard   

  end  % igra_ef 

else
%  n_ve1=n_ve;
%  n_pa1=n_pa; 
end

 Delta_eps=(ex-ey)/2;
 %' igrapp', keyboard 
% Delta_eps=(ex-ey)/4;

if ifp==-10
 ' igrapp', keyboard
end
 if igr_app==1 & igraef_new<=1

  fi6=find(shailoop~=6);
  shailoopu=shailoop(fi6);
  sha=0;
  rad(ind_vet,:)=0;
  aitot(ind_vet,:)=0;
  nv(ind_vet,:)=0;
  if iLP==0
   nv(ind_vet,1)=sqrt(n_me);
   nitot(ind_vet-1,1)=sqrt(n_me);
%   anyf(ind_vet,:)=1;
%    nv(ind_vet,1)=sqrt(ex);
%    nv(ind_vet,1)=sqrt(ey);
  else
   if icambio_ref==1
    nv(ind_vet,1)=sqrt(ex);
   else
    nv(ind_vet,1)=sqrt(ey);
   end
  end
 elseif igr_app==2
  if iLP==0
   nv(ind_vet,1)=sqrt(n_me);
%'wuo', keyboard
  else
   if icambio_ref==1
    nv(ind_vet,1)=sqrt(ex);
   else
    nv(ind_vet,1)=sqrt(ey);
   end
  end
  sha=1;
  shag=sha;
  lxi=1;
  rgra=P.Ry;
  adis=rgra*kcav;
  rad(ind_vet,:)=rgra;
  aitot(ind_vet,:)=rgra;
%'wuo  las' , keyboard
%  nv(ind_vet,2)=sqrt(max([er1 er2]));
  nv(ind_vet,2)=sqrt(er2);

%    nv(ind_vet,1)=sqrt(ex);
%    nv(ind_vet,1)=sqrt(ey);

 end
 
 shavet(ind_vet)=sha;

%'set shavet', keyboard
%'set shavet', 

 nitn=nv(put,:).';
 nib =nv(pub,:).';
 shav.t=shavet(put,:);
 shav.b=shavet(pub,:);
 if iLP==0
  anyret(ind_vet)=2*Delta_eps;
 end
 anyr.t=anyret(put,:);
 anyr.b=anyret(pub,:);
 aitn=rad(put,:).';
 aib=rad(pub,:).';

% bi=aral.y.t';
% pai=aral.p.t';
if ifp==-10
'prima di ireturn =1 ', keyboard
end
 if igraef_new<2
  ireturn=1;
  inu_int=1
  if inu_int==1
  rado=rand;
  ran=num2str(rado*1e20);
  fipun=find(ran~='.');
  loc=['loc',ran(fipun(1:4))];
  
%  eval(['save ',loc])

save n1d
%'in shagrat ', keyboard
%aysa=ay;
%fi2d  
%ay=aysa;

if ifp==-10
'in shagrat tolto ricalcolo campo ', keyboard
end
 %     lambdan=mean([la_pa la_ve]);
      
%      'in shagrat dopo lambdan ', keyboard

      
%  eval(['load ',loc])
%  eval(['!del ',loc,'.mat'])
%  return
  

inochaf=1;
if inochaf==0
  fre_ca=linspace(fre_camp(1),fre_camp(end),fix(Ndisp/2));
  fre_camp=[fre_ca fre_ca+De_fre];

%  fre_ca=linspace(fre_camp(1),fre_camp(end),Ndisp);
%  fre_camp=[fre_ca+De_fre];  
%  fre_camp=[fre_ca fre_ca+De_fre]-De_fre/2;
%   kcav=2*pi/(lambda-De_lam/2)*rr;
  
 Ndisp=-1;
 
losalto=0
if losalto==0
 if mm==0 & mmvet(1)~=0
   fre_camp=fre_camp+Dlam_mod(4)/lambda/1000;
  end
end
%fre_camp=linspace(Dlam_mod(1),Dlam_mod(2),11)/lambda/1000;
if ifp==-10
 figure, plot(fre_camp*1000*lambda,'.'), keyboard
end
  freq=fre_camp(1);
'sistemo inter', keyboard
end  %inochaf  
  
%'sistemo inter', keyboard  

%'sistemo inter', keyboard
  end
  if ifp~=-4
   ' SHAGRAT: fine calcolo ret efficace', keyboard
   close all
  end
  return
 end
end


ilogra=0
if ilogra==1
 load retf
 return
 % save retf rv AB A B
end 
if ~exist('Rma')
Rma=1
end

if Rext>Rma
 Rvero=Rext;
 shiff=1;
else
 shiff=P.shif;
 Rvero=Rma;
end

if shiff~0
 shiret=-D_stri/2;
else
 shiret=0;
end

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


   switch Shap
    case 'circle'

     Rma=Y;
     cas=0;
    case {'square','rectangle'}
     cas=1;
     Cmae=(ns-1)/2*D_stri+d_stri/2;
     if orgra==0
      Rma=sqrt(Cmae^2+X^2);
     else
      Rma=sqrt(Cmae^2+Y^2);
     end
   end


%npr=201*ceil(Rma/D_stri);
%rv=linspace(0,Rma*1.01,npr);

npr=201*fix(Rvero/D_stri);
rv=linspace(0,Rvero,npr);
if Rext>0
 fireti=find(rv<=Rma);
 rv_ret=rv(fireti);
else
 rv_ret=rv;
end


if ifp>-3 | ifp==-10 | ifp>=0

  figure

iold=0;
 if iold==1
%  yv0=[-d_stri -d_stri d_stri d_stri]/2;
  yv0=[-d_stri*ones(1,100) d_stri*ones(1,100)]/2;
  for k=1:ns
   Ym=(k-1-fix(ns/2))*D_stri-shiret;
   yvd=yv0+Ym;
    if cas==0
     dif=Xlo^2-(Ym)^2;
    elseif cas==1
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

 else


%npr=201*fix(Rma/D_stri);
%rv=linspace(0,Rma,npr);


  yv0=d_stri/2*[-1 1];
  npf=50;
  rtot=[];
  rtoti=[];
  for k=1:ns
   Ym=(k-1-fix(ns/2))*D_stri-shiret;
   yvd=yv0+Ym;
   aas=yvd/Rma;
   fias=find(abs(aas)<=1);
   fiv=asin(aas(fias));
   if length(fiv)~=0
    if length(fiv)==1
     if fiv<0
      five=linspace(pi-fiv,2*pi+fiv,2*npf);
     else
      five=linspace(fiv,pi-fiv,2*npf);
     end
     five=five+orgra;
     rfiu=Rma*exp(j*(five));
     fiex=Rma*cos(five([1 end]));
     rfiix=linspace(fiex(1),fiex(2),2*npf);
     rfiiy=ones(size(rfiix))*Rma*sin(five(1));
     rfii=rfiix+j*rfiiy;
     rtoti=[rfii rfiu];
     rtot=[rtot; rtoti];
    else
     five1=linspace(fiv(1),fiv(2),npf)+orgra;
     five2=linspace(pi-fiv(2),pi-fiv(1),npf)+orgra;
     rfiu1=Rma*exp(j*five1);
     rfiu2=Rma*exp(j*five2);
     fiex=Rma*cos([five1(end) five2(1)]);
     rfiix=linspace(fiex(1),fiex(2),npf);
     rfiiy=ones(size(rfiix))*Rma*sin(five1(end));
     rfii1=rfiix+j*rfiiy;

     fiex=Rma*cos([five2(end) five1(1)]);
     rfiix=linspace(fiex(1),fiex(2),npf);
     rfiiy=ones(size(rfiix))*Rma*sin(five1(1));
     rfii2=rfiix+j*rfiiy;
     rtoti=[rfiu1 rfii1 rfiu2 rfii2];
     rtot=[rtot; rtoti];
    end
   end
  end
    if r_add>0
     fiad=find(abs(rtot)<r_add);
     rtot(fiad)=NaN;
    end

    plot3(real(rtot.'),imag(rtot.'),ones(size(rtot.')),'w'), view(2), hold on

 end

    if r_add>0
     fi=linspace(0,2*pi,155);
     cirv=r_add*exp(j*fi);
     plot3(real(cirv),imag(cirv),ones(size(cirv)),'w'), view(2)
    end

%  npr=1001;

  hold on,
   switch Shap
    case 'circle'
     polar(linspace(0,2*pi,100),max(rv)*ones(1,100),'r')
     if Rext>0
      polar(linspace(0,2*pi,100),Rma*ones(1,100),'r')
      title(' the areas within white and red lines is etched away ')
     end
    case {'square','rectangle'}
     if orgra==0
      xplo=[X X -X -X X];
      yplo=[-1 1 1 -1 -1]*Cmae;
     else
      xplo=[1 1 -1 -1 1]*Cmae;
      yplo=[-Y Y Y -Y -Y];
     end
     plot(xplo,yplo,'r')
   end
  axis equal
  drawnow
  if ifp>-3
   pausak
  end
%'cirv in sha_grat'  , keyboard
end

%if ifp~=-4

%end
%muv=[0:2:2*(nubesu+1)];
PAR_sha.nubesu=nubesu;
PAR_sha.r_add=r_add;
PAR_sha.Rma=Rma;
PAR_sha.rv_ret=rv_ret;
PAR_sha.rv=rv;
PAR_sha.shiret=shiret;
PAR_sha.cas=cas;
PAR_sha.d_stri=d_stri;
PAR_sha.D_stri=D_stri;
PAR_sha.X=X;
PAR_sha.Y=Y;
PAR_sha.Ylo=Ylo;
PAR_sha.orgra=orgra;
PAR_sha.kcav0=kcav0;
PAR_sha.pimu=pimu;
PAR_sha.meun=meun;
PAR_sha.mbv=mbv;
PAR_sha.iclo=iclo;


%' verifica dove sono prima function', keyboard
[A,B,AB]=sha_grat_fun(PAR_sha,ifp);
r=rv*kcav0;
rv=r;
%' verifica finito', keyboard
return

    if r_add>0
     firv=find(rv_ret>r_add & rv_ret<=Rma);
     firv0=1:firv(1)-1;
    else
     firv=find(rv_ret>0 & rv_ret<=Rma);
     firv0=1;
    end

%AB=ones(length(rv),length(muv))*NaN;
AB=zeros(length(rv),length(muv));
AB(:,1)=pi;
im=0;
for mu=muv

im=im+1;
%mu, keyboard
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


   if cas==0
   elseif cas==1
   end

  F_sh(firv0)=fat1;
%  firv=100;
  for k=firv
%   if k==1000
%   etis3=etime(clock,tc1)
%   pausak
%   end
%k
%tic
   r=rv(k);
   if cas==0
    a0=d_stri/(2*r);
   elseif cas==1
    H1=sqrt(X^2+(d_stri/2)^2);
    if r<=X
     a0=d_stri/(2*r);
    elseif r>X & r<=H1
     H=sqrt(r^2-Ylo^2);
     a0=H/r;
    elseif r>H1
     a0=0;
    end
   end
   nstr=fix((r+d_stri/2-shiret)/D_stri);
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
    if cas==0
     a1=(kj*D_stri+d_stri/2+shiret)/r;
    elseif cas==1
     Y2=kj*D_stri+d_stri/2+shiret;
     H2=sqrt(X^2+(Y2)^2);
     if r<=H2
      a1=(kj*D_stri+d_stri/2+shiret)/r;
     elseif r>H1 & r<=H2
 %     H=sqrt(r^2-H1^2);
      H=sqrt(r^2-Ylo^2);
      a1=H/r;
     elseif r>H2
      a1=0;
     end
    end
    if a1<=1
     as1=asin(a1);
    else
     as1=fat;
    end

    if cas==0
     a2=(kj*D_stri-d_stri/2+shiret)/r;
    elseif cas==1
     Y1=kj*D_stri-d_stri/2+shiret;
     Y2=kj*D_stri+d_stri/2+shiret;
     H1=sqrt(X^2+(Y1)^2);
     H2=sqrt(X^2+(Y2)^2);
     if r<=H1
      a2=(kj*D_stri-d_stri/2+shiret)/r;
     elseif r>H1 & r<=H2
%      H=sqrt(r^2-H1^2);
      H=sqrt(r^2-Ylo^2);
      a2=H/r;
     elseif r>H2
      a2=0;
     end
    end
    as2=asin(a2);
    som=som+(as1-as2);
   end

   F_sh(k)=som+som0;
   if imag(F_sh(k))~=0
%    pausak
   end
%  toc
%  pausak
  end  %r
  if ifp>0
%   ver=4*(F_sh*rv')*diff(rv(1:2))/(pi*max(rv)^2)/(d_stri/D_stri)
%   pausak
  end
  F=F_sh;
 else

  fat=0;
  F_sh(1)=fat;
  if r_add>0
   F_sh(firv0)=pi/4; %ERRORE
   F_sh(firv0)=0;
  end
  for k=firv
   r=rv(k);
   if cas==0
    a0=d_stri/(2*r);
   elseif cas==1
    H1=sqrt(X^2+(d_stri/2)^2);
    if r<=X
     a0=d_stri/(2*r);
    elseif r>X & r<=H1
%     H=sqrt(r^2-(d_stri/2)^2);
     H=sqrt(r^2-Ylo^2);
     a0=H/r;
    elseif r>H1
     a0=0;
    end
   end
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
    if cas==0
     a1=(kj*D_stri+d_stri/2+shiret)/r;
    elseif cas==1
     Y2=kj*D_stri+d_stri/2+shiret;
     H2=sqrt(X^2+(Y2)^2);
     if r<=H2
      a1=(kj*D_stri+d_stri/2+shiret)/r;
     elseif r>H1 & r<=H2
%      H=sqrt(r^2-H1^2);
      H=sqrt(r^2-Ylo^2);
      a1=H/r;
     elseif r>H2
      a1=0;
     end
    end

    if a1<=1
     as1=asin(a1);
    else
     as1=fat;
    end
    if cas==0
     a2=(kj*D_stri-d_stri/2+shiret)/r;
    elseif cas==1
     Y1=kj*D_stri-d_stri/2+shiret;
     Y2=kj*D_stri+d_stri/2+shiret;
     H1=sqrt(X^2+(Y1)^2);
     H2=sqrt(X^2+(Y2)^2);
     if r<=H1
      a2=(kj*D_stri-d_stri/2+shiret)/r;
     elseif r>H1 & r<=H2
%      H=sqrt(r^2-H1^2);
      H=sqrt(r^2-Ylo^2);
      a2=H/r;
     elseif r>H2
      a2=0;
     end
    end
    as2=asin(a2);
    som=som+(sin(mu*(as1+orgra))-sin(mu*(as2+orgra)))/mu;
   end
   F_sh(k)=som+som0;
  end  %r
%  F_sh=F_sh;
  F=F_sh;
 end
 AB(1:firv(end),im)=2*F_sh';
 clear F_sh
%etis2=etime(clock,tc1)
 if ifp>=2
  im
 end
%pausak

% figure, plot(rv,F_sh,'.g')
% pausak
 if ifp>=2
  if exist('F')
   figure, plot(rv,F,'.g')
   pausak
  end
 end
 clear F
end %im

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
' dopo primo passo', keyboard
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
%   [mu nu], pausak
  end
 end

 if ifp>=2
  ' fine ahgrat', keyboard
 end
%a=asav;
disp(' fine sha_grat '), keyboard
%figure, plot(rv/a*avero,AB), keyboard
etis=etime(clock,tc1)
%if ifp~=-4
' fine sha_grat'
   if iclo==1
    close
   end
% keyboard
%end
