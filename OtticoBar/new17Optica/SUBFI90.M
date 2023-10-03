if ~exist('kcut')
 kcut=0;
end
%' usat', keyboard

 shi=length(Aue1)/2;
 Kfi=Kad;
 puki=puk(1):puk(2);
 ad=Aue1(puki);
 KKv0=KKtdu(puki);
 if ialt==1
  Aued=[0; ad];
  KKv0=[0; KKv0];
 else
  Aued=ad;
 end
%' usat 0', keyboard
 Aspl=spline(KKv0,Aued,Kfi);
%' usat 1', keyboard
 if kcut>0
  fik=find(Kfi>kcut);
  Aspl(fik)=0;
 end
% figure, plot(KKv0,abs(Aued),Kfi,abs(Aspl),'o'), pausak
% figure, plot(KKv0,real(Aued),Kfi,real(Aspl),'o'), pausak
% figure, plot(KKv0,imag(Aued),Kfi,imag(Aspl),'o'), pausak
 Aue21=[Aue21; Aspl];

%' ubom', keyboard
if iLP==0
 ad=Aue1(puki+shi);
% KKv0=KKtdu(puki+shi);
 if ialt==1
  Aued=[0; ad];
 else
  Aued=ad;
 end

 Aspl=spline(KKv0,Aued,Kfi);
 if kcut>0
  fik=find(Kfi>kcut);
  Aspl(fik)=0;
 end
 Aue22=[Aue22; Aspl];

% ez

 ad=Auez1(puki);
 if ialt==1
  Aued=[0; ad];
 else
  Aued=ad;
 end
 Aspl=spline(KKv0,Aued,Kfi);
 if kcut>0
  fik=find(Kfi>kcut);
  Aspl(fik)=0;
 end
 Auez21=[Auez21; Aspl];
% Auez21=[Auez21; spline(KKv0,Aued,Kfi)];

 ad=Auez1(puki+shi);
 if ialt==1
  Aued=[0; ad];
 else
  Aued=ad;
 end
 Aspl=spline(KKv0,Aued,Kfi);
 if kcut>0
  fik=find(Kfi>kcut);
  Aspl(fik)=0;
 end
 Auez22=[Auez22; Aspl];
% Auez22=[Auez22; spline(KKv0,Aued,Kfi)];

% hz

 ad=Auhz1(puki);
 if ialt==1
  Aued=[0; ad];
 else
  Aued=ad;
 end
 Aspl=spline(KKv0,Aued,Kfi);
 if kcut>0
  fik=find(Kfi>kcut);
  Aspl(fik)=0;
 end
 Auhz21=[Auhz21; Aspl];
% Auhz21=[Auhz21; spline(KKv0,Aued,Kfi)];

 ad=Auhz1(puki+shi);
 if ialt==1
  Aued=[0; ad];
 else
  Aued=ad;
 end
 Aspl=spline(KKv0,Aued,Kfi);
 if kcut>0
  fik=find(Kfi>kcut);
  Aspl(fik)=0;
 end
 Auhz22=[Auhz22; Aspl];
% Auhz22=[Auhz22; spline(KKv0,Aued,Kfi)];

% h

 ad=Auh1(puki);
 if ialt==1
  Aued=[0; ad];
 else
  Aued=ad;
 end
 Aspl=spline(KKv0,Aued,Kfi);
 if kcut>0
  fik=find(Kfi>kcut);
  Aspl(fik)=0;
 end
 Auh21=[Auh21; Aspl];
% Auh21=[Auh21; spline(KKv0,Aued,Kfi)];

 ad=Auh1(puki+shi);
 if ialt==1
  Aued=[0; ad];
 else
  Aued=ad;
 end
 Aspl=spline(KKv0,Aued,Kfi);
 if kcut>0
  fik=find(Kfi>kcut);
  Aspl(fik)=0;
 end
 Auh22=[Auh22; Aspl];
% Auh22=[Auh22; spline(KKv0,Aued,Kfi)];

  bespd=besselj(ibp,Kad*xvero*kref);
  beszd=besselj(ibz,Kad*xvero*kref);
  Mvegmd=repmat(Mvegm0(puk(1),if90),Np,1);
  Mvefpd=repmat(Mvefp0(puk(1),if90),Np,1);
  Mvegpd=repmat(Mvegp0(puk(1),if90),Np,1);
  Mvezd=repmat(Mvez0(puk(1),if90),Np,1);
  Mvhzd=repmat(Mvhz0(puk(1),if90),Np,1);
%  [kin ialt], pausak

  besp=[besp; bespd];
  besz=[besz; beszd];
  Mvefp=[Mvefp; Mvefpd];
  Mvegm=[Mvegm; Mvegmd];
  Mvegp=[Mvegp; Mvegpd];
  Mvez =[Mvez;  Mvezd ];
  Mvhz =[Mvhz;  Mvhzd ];

end

  besmd=besselj(ibm,Kad*xvero*kref);
  Mvefmd=repmat(Mvefm0(puk(1),if90),Np,1);
  besm=[besm; besmd];
  Mvefm=[Mvefm; Mvefmd];

%' las', keyboard
