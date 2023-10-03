%' ferma', keyboard
iome=0;
if iome==1
 fsr_sub
else 

global fiload Ppol uL
%eval([' save ',fiload])
% fiload='sol0';
 fisha=find(shavet(:,1)==4);
 if length(fisha)==0
 fisha=find(shavet(:,1)<4);
 end
 aimp=radii.a(fisha);
 Rmai=max(aimp);
 
 R0v=Rmai*linspace(0,1.2,11); 
 R0v=0; 
%'fiload;\\/\',  keyboard
% R0v=Rmai*0.8; 
 kvv=0;
% kvv=Ps.KKmax;
% kvv(1)=.00;
 iraffina=1;
% iraffina=0;

ive=0;
clear Ppol

if ive==1
  R=R0v; 
%' nuomo posto', keyboard
 for ke=1:length(kvv)
 kv=kvv(ke);
 [lambda2,uL,nto,dto,gth2]=sol1_fun(R0v,fiload,iraffina,kv);


Ppol.uL=uL;
%' nuomo posto', keyboard
% fi=find(gth2>0);
% R=R0v(fi);
% L=lambda2(fi);
% G=gth2(fi);


  L(:,ke)=lambda2';
  G(:,ke)=gth2';
end

else  %ive
 iref=iraffina;
 R0=0;
 kkee=1;
 ks=kkee;

%' sol0', keyboard 
 sol0_sub

%' dopo sol0', keyboard 

 if kkee>=2
  ifp=-4;
 end
  if kkee==1
   nto1=nto;
   dto1=dto;
  end
 if iany~=3 
  lambda2(1)=lambda;
  gth2(1)=gth;
  if ~exist('lambda1')
   lambda1=lambda;
   gth1=gth;
  end 
  lambda2(2)=lambda1;
  gth2(2)=gth1;
  Ppol.uL=uL;
  L(:,ks)=lambda2';
  G(:,ks)=gth2';
 else
  Ppol.uL=uL;

  if ipolar==-1
   lambda2(1)=lambda;
   gth2(1)=gth;  
   L(:,1)=lambda;
   G(:,1)=gth; 
  else
   lambda2(2)=lambda;
   gth2(2)=gth;
   L(:,2)=lambda;
   G(:,2)=gth;    
  end
 end
%' nuovo ', keyboard
end  %ive



Le=max(L);
L_sup=max(max(L));
L_inf=min(min(L));
Dlam=(L_sup-L_inf)*1000;
lame=L_sup;

lambda=lambda2(1);
gth1d=gth2(1);
if i1D==1
 return
end

iBel=1;
if iBel==0
Gm=mean(mean(G));
fi=find(G<Gm);
L_sup1=max(max(L));
L_inf1=min(min(L));  
Lm=mean(L(fi));
Dlam=(Lm-L_inf)*1000;
lame=Lm;
Dlam_mod(1)=-0.2;
end

%lambda=lame;

%if Dlam<0.5
% Dlam_mod(2)=0.5;
%else
% Dlam_mod(2)=Dlam;
%end
%Dlam_mod(2)=Dlam_mod(2)+Dlam;

%'ferma dopo', keyboard

%'dopo BEL', keyboard
isl=1;

if isl==0
 lambda=max(lambda2);
 ls=sort(lambda2);
 dlambda=diff(ls([1 end]))*1000;
 dlambda=0
     if dlambda>Dlam_mod(2)
      Dlam_mod(2)=dlambda;
      Dlam_mod(1)=0;
     end 
     if Dlam_mod(5)>0
      NP=fix(diff(Dlam_mod(1:2))/Dlam_mod(5));
%      ' NP ', keyboard
      Dlam_mod(3)=NP;
     end 
end
 if ifp==-10 & length(R0v)>1
  R=R0v;
  O=ones(size(R0v));
  figure
  subplot(211), plot(R,L*1000,'.-',R,O*L_inf*1000,'r',R,O*L_sup*1000,'g')
  xlabel([' radial coord. (micron)'])
  subplot(212), plot(R,G),
  pausak
 end
 
end 

%' fine 1D ', keyboard