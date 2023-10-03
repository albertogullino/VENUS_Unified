close all
addpath E:\Dati\mvcsel\new10Optica
rmpath E:\Dati\mvcsel\new10Vortex
colordef black

MGiri=Par.Ngiri;
o=ones(size(KK));
ot=[];
for km=1:numodi
 ot=[ot; o*(km-1)];
end

over=[ot; ot; ot; ot];

%figure, plot([1:length(over)],over), 
%hold on, plot(Pus0,over(Pus0),'r.'), pausak

Ps.isav_Az=0

nma=1;
isolut=0;




%icampi=0
%ad
%close all
isk=1
if isk==0
if exist('R0v')
    R=R0v;
  O=ones(size(R0v));
Gm=mean(mean(G));
fi=find(G<Gm);
L_sup1=max(max(L));
L_inf1=min(min(L));  
Lm=mean(L(fi));
  figure
  subplot(211), plot(R,L*1000,'.-',R,O*L_inf*1000,'r',R,O*L_sup*1000,'g')
  hold on, plot(R,O*L_inf1*1000,'r--',R,O*L_sup1*1000,'g--',R,O*Lm*1000,'c-o')
  xlabel([' radial coord. (micron)'])
  subplot(212), plot(R,G,R,O*Gm,'r')
  pausak
end
end
clear Gsov Fsov
isavetutto=0;
%ipolar=2;
nmasce=1;
%'postg', keyboard
ipostg=1;        %=0 salva
ipost=0;        %=0 ricalcola base
%ipost=1;        %=0 ricalcola base
%isav_Az=0;        %=0 calcolo campi z
%Ps.isav_Az=isav_Az;        %=0 calcolo campi z
%ipost=0;        %=0 ricalcola base
iade=0;    % fa come in programma
%clear Im
ifp=-10;
%ifp=-4;
%icam_fr=input(' campi per ogni frequenza ? [0/1] ');
if ~exist('icam_fr')
 icam_fr=0;
end
%ionly1=input(' solo una soluzione ? [0/1] ');
ionly1=1;
if length(ionly1)==0
 ionly1=0;
else

% imod=input(' modo 1 o 2 ? [1/2] ');
% imod=1;
 ICON=0;
% ICON=0;
  ipolar=pvet;
end
ipolarv=length(ipolar);
%clear global
%load
%ipolar=-1;
%  ifp=-10;
%  if exist('Pus0')
%   Pus=Pus0;
%  end


Over=over(Pus0);

Mazv{1}=[100];
mmm=[0:20];
mmm0=MGiri+[-1 1];
fi1=find(mmm~=mmm0(1) & mmm~=mmm0(2));
mmmd=mmm(fi1);
Mazv{2}=mmm0;
Mazv{3}=mmmd;
MaxArm=3;

%ha=figure;
%set(ha,'pos',[624          10        1078         974]);

 ipost=0;

 dissfun='diss_Vorarm';
 
  ifps=ifp; if ifp>=-2; ifp=1; end
     pola=pvet(1);
     Azvet=Azvet1;
     Azvetf=Azvetf1;
     Gvet=Gvet1;
   alvet=alvet1;
   
  if  length(desi)==0
   desi=input(' plot a sinistra [1] o destra [0] = ');
   if length(desi)==0
    desi=0;
   end
  end
  enlar=0;
  nsub2=1;
  enlar=300;
  nsub2=2;
    ved=version;
    ve=ved(1);
   if desi==1
    sinc=+1;
    pograp=[ 40 70 450 600];
    pogram=[ 40 30 450 600];

    if str2num(ve(1))==7
     pograp=[33   200   755   772];
     pogram=[33   104   755   772];
     pograp=[20   200   500   550];
     pogram=[20   100   500   550];
    end 
    else
    sinc=-1;
    pograp=[ 550 70 450 600];
    pogram=[ 550 30 450 600];

    if str2num(ve(1))==7
     pograp=[550   0   500   600];
     pogram=[550   50   500   600];
    end 
   end
   


for imazver=1:MaxArm

imazver
Maz=Mazv{imazver};

Ofin=zeros(size(Over));
fiz=[];
for kkk=1:length(Maz)
 fiz=[fiz; find(Over==Maz(kkk))]; 
end

if length(fiz)==0
 Ofin=ones(size(Over));
else 
 Ofin(fiz)=1;
end







  imod=0;
  isalva=1;

  camv_arm
 %if imazver>1  
 % close(4:5)
 %end 
end


ha1=figure;
set(ha1,'pos',[593         612        1313         367]);
etot=sum(effintv(2:end));
s3=size(EF,3);
for imazver=1:s3
             figure(ha1)
             subplot(1,3,imazver)
             surf(X,Y,EF(:,:,imazver)/effma), colorbar
%             surf(X,Y,real(EF(:,:,imazver))/sqrt(effma)), colorbar
                         shading('interp'), view(0,90),
	                 axis square, axis equal, grid,
	                 axis([-1 1 -1 1]*axlim/2),
	                 if imazver==1
	                  tit=[' Totale; Integrale =1 '];
	                 elseif imazver==2
	                  tit=[' Purezza Integrale =',num2str(effintv(imazver)/etot,3)];
	                 elseif imazver==3
	                  tit=[' Contibuti spuri =',num2str(effintv(imazver)/etot,3)];	                  
	                 end
	                 title(tit)
            
            
end            