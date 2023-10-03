%Lvtv=[d*1e6; Lvt];
%nvtv=[rr; nvt];


   [Gu0,du,Trasup0]=gaemms(0,freq,lambda,Lplit,nplit,Lplt,nplt,Nsplt,...
              rfu,rr,iLP1,Lplut,nplut,rqw);
  [Gu1,du,Trasup1]=gaemms(0,freq,lambda,Lplit,real(nplit),Lplt,real(nplt),...
              Nsplt,rfu,rr,iLP1,Lplut,real(nplut),rqw);


% Lplit1=flipud(Lplit);
% Lplut1=flipud(Lplut);
% Lplt1=flipud(Lplt);
% nplit1=flipud(nplit);
% nplut1=flipud(nplut);
% nplt1=flipud(nplt);

% [Lplit nplit; 0 0; Lplt nplt; 0 0; Lplut nplut]

%  [Gu0,du,Trasup0]=gaemms(0,freq,lambda,Lplit,real(nplit),Lplt,...
%              real(nplt),Nsplt,rfu,rr,iLP1,Lplut,real(nplut),rqw);
  [Gu1,du,Trasup1]=gaemms(0,freq,lambda,Lplit,(nplit),Lplt,...
              (nplt),Nsplt,rfu,rr,iLP1,Lplut,(nplut),rqw);
%  [abs(Trasup0)^2 abs(Gu0)^2 1-abs(Gu0)^2-abs(Trasup0)^2]
%  [abs(Trasup1)^2 abs(Gu1)^2 1-abs(Gu1)^2-abs(Trasup1)^2]


%  [Gb0,du,Trainf0]=gaemms(0,freq,lambda,Lplib,real(nplib),Lplb,real(nplb),...
%              Nsplb,real(rfd),rr,iLP1,Lplub,real(nplub),rqw);
  [Gb1,du,Trainf1]=gaemms(0,freq,lambda,Lplib,nplib,Lplb,nplb,Nsplb,...
              rfd,rr,iLP1,Lplub,nplub,rqw);
%  [abs(Trainf0)^2 abs(Gb0)^2 1-abs(Gb0)^2-abs(Trainf0)^2]
%  [abs(Trainf1)^2 abs(Gb1)^2 1-abs(Gb1)^2-abs(Trainf1)^2]

iff=1;  %per calcolare il campo in z

iff=0;  %per calcolare il campo in z

%'gaz7 prima', keyboard
if ifp==-10
%'gaz7 prima', keyboard
end
   [gazk,Iku,Ikb,fak,Tru,Trb,Gu,Gb,Le,Lcav,ztot,Ez,Hz,indz,nmean,Perd]=...
   gaz7(fiQ,fiCav,L_i,n_i,rr,rfd,rfu,lambda,freq,0,iLP,ifp,iff);

%  [abs(Tru)^2 abs(Gu)^2 1-abs(Gu)^2-abs(Tru)^2]
%  [abs(Trb)^2 abs(Gb)^2 1-abs(Gb)^2-abs(Trb)^2]

%  [gazk,Iku,Ikb,fak,Trur,Trbr,Gur,Gbr,Le,Lcav,ztot,Ez,Hz,indz,nmean]=...
%  gaz7(fiQ,fiCav,L_i,real(n_i),rr,real(rfd),rfu,lambda,freq,0,iLP,ifp,iff);
%  [abs(Trur)^2 abs(Gur)^2 1-abs(Gur)^2-abs(Trur)^2]
%  [abs(Trbr)^2 abs(Gbr)^2 1-abs(Gbr)^2-abs(Trbr)^2]

Lemi=abs(Le);
Lef=abs(Le*1e-6);
iofa=1;
if iofa==1
 fatqw=fak;
 confz=abs(gazk/(fatqw*NQW));
 confztot=abs(gazk);
else
 confz=confzv(1);
 confztot=confzv(2);
end

fTras=abs(Tru)^2;
fTrasinf=abs(Trb)^2;
fPsup=1-abs(Gu)^2-fTras;
fPinf=1-abs(Gb)^2-fTrasinf;


%disp(' T_ef '), keyboard
%disp(' T_ef '), keyboard
%disp(' T_ef '), keyboard

if iff==1
% figure, plot(ztot,abs(Ez).^2*3,ztot,real(indz))
% figure, plot(ztot,abs(Ez).^2*3,ztot,real(indz))
 figure, plot(ztot,abs(Ez).^2*3,ztot,abs(Hz).^2*3,ztot,real(indz))
 if ifp>-4
  pausak
 end
end
%disp(' tl_ef '), keyboard