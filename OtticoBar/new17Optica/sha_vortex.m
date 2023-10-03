%ifp=10
%if istrumix==1
% aS=a;
% bS=b;
%
% a=aloc;
% b=bloc;
%end
%igint=0;
igint=1;
' sha_vortex', keyboard
iplo=1;
if ifp==-4
 iplo=0;
end 

clear A B C D CD AB

mazim=Par.mazim;
SogDen=Par.SogDen;

%'ot', keyboard


an1=acos(SogDen);
an2=acos(-SogDen);

pea=2*pi/mazim;

anv0=[an1 an2]/mazim;
anv=[anv0 pea-fliplr(anv0)];
sogv=[SogDen -SogDen -SogDen SogDen];


an_vet=[];
so_vet=[];
for k=1:mazim
 an_vet=[an_vet anv+(k-1)*pea];
 so_vet=[so_vet sogv];
end

an_veti=an_vet([end 2:2:end-2]);
an_vetu=an_vet(1:2:end-1);

[an_veti' an_vetu']

if iplo==1
firad=linspace(0,2*pi,501);
figure, plot(firad,cos(mazim*firad),an_veti,cos(mazim*an_veti),'go',an_vetu,cos(mazim*an_vetu),'ro'), 
hold on
plot(firad,ones(size(firad))*SogDen,'w--',firad,-ones(size(firad))*SogDen,'w--')
pausak
%figure, plot(firad,sin(mazim*firad),an_veti,sin(mazim*an_veti),'go',an_vetu,sin(mazim*an_vetu),'ro'), 
%hold on
%plot(firad,ones(size(firad))*SogDen,'w--',firad,-ones(size(firad))*SogDen,'w--')
%pausak
end
'ferma vor'
keyboard
  aa=adis;
  ab=bdis;
  
  if igint==1
   ng=101;
   [rv,wi]=gauleg(aa,ab,ng);
  else
   rv=linspace(aa,ab,ng);
  end
  r=rv/(kcav);  
  sgimp=ones(size(r));

muv=[0:pasnu:2*nubesu];    
AB=ones(2,length(muv))*NaN;
CD=ones(2,length(muv))*NaN;
im=0;
ab0d=an_vetu-an_veti;  
ab0=ab0d(2)*mazim*2;

for mu=muv
im=im+1;
 if mu==0
  AB(1:length(rv),im)=ab0/2;
  CD(1:length(rv),im)=0;
 else
  AB(1:length(rv),im)=1/(2*mu)*sum(sin(mu*an_vetu)-sin(mu*an_veti));
  CD(1:length(rv),im)=-1/(2*mu)*sum(cos(mu*an_vetu)-cos(mu*an_veti));
 end
 %(sin(mu*an_vetu)-sin(mu*an_veti))', pausak
end

fi=find(abs(AB)<1e-10);
AB(fi)=0;

fi=find(abs(CD)<1e-10);
CD(fi)=0;
%'qui',  keyboard

rodis=r;
if iplo>2
 if length(muv)>6
  figure, plot(rodis,AB(:,1:6)), hold on, plot(rodis,AB(:,7:length(muv)),'.-'),
 else
  figure, plot(rodis,AB),
 end
 title(' AB: y0, m2, c4, r6, g8, b10, y12, m14, c16, r18, g20, b22')
 pausak
 if length(muv)>6
  figure, plot(rodis,CD(:,1:6)), hold on, plot(rodis,CD(:,7:length(muv)),'.-'),
 else
  figure, plot(rodis,CD),
 end
 title(' CD: y0, m2, c4, r6, g8, b10, y12, m14, c16, r18, g20, b22')
 pausak
end

  if iplo>10
   figure
  end
 for imu=pimu
  jmu=imu-meun;
  mu=mbv(imu);
  for inu=pimu
   jnu=inu-meun;
   nu=mbv(inu);
%   if (nu+mu)/2-fix((nu+mu)/2)==0
    dmn=abs(mu-nu);
    sdmn=sign(mu-nu);
%    if sdmn==0
%     keyboard
%    end
    fim=find(dmn==muv);
    if length(fim)==1
     mfatd=AB(:,fim);
     A(:,jmu,jnu)=AB(:,fim);
     if ipolar==0
      C(:,jmu,jnu)=sdmn*CD(:,fim);
     end
    else
     disp('errore A mu in sha_oxi ')
     pausak
    end
    dmn=abs(mu+nu);
    fim=find(dmn==muv);
    if length(fim)==1
     mfats=AB(:,fim);
     B(:,jmu,jnu)=AB(:,fim);
     if ipolar==0
      D(:,jmu,jnu)=CD(:,fim);
     end
    else
     disp('errore B mu in sha_oxi ')
     pausak
    end
   end  %if
%  end
 end
 if iplo>=1
  figure, plot(rodis,AB), hold on, plot(rodis,CD,'--')
  title(' AB continuous, CD dashed')
  keyboard
 end  
  
  
return
muv=[0:pasnu:2*nubesu];
%muv=[0:2:8*nubesu];
%muv=[0:2:18];
AB=ones(2,length(muv))*NaN;
im=0;
ab0=mazim*(an_vetu(2)-an_veti(2));
for mu=muv
im=im+1;
 if mu==0
  AB(1:length(rv),im)=ab0;
 else
  AB(1:length(rv),im)=1/(2*mu)*sum(sin(mu*an_vetu)-sin(mu*an_veti));
 end
 %(sin(mu*an_vetu)-sin(mu*an_veti))', pausak
end

fi=find(abs(AB)<1e-10);
AB(fi)=0;

%'qui',  keyboard

clear A B
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
     disp('errore A mu in sha_oxi ')
     pausak
    end
    dmn=abs(mu+nu);
    fim=find(dmn==muv);
    if length(fim)==1
     mfats=AB(:,fim);
     B(:,jmu,jnu)=AB(:,fim);
    else
     disp('errore B mu in sha_oxi ')
     pausak
    end
   end  %if
  end
 end

if iplo==1
% figure, plot(r,A(:,1:6)), hold on, plot(r,A(:,7:length(muv)),'.-'),
% title(' y0, m2, c4, r6, g8, b10, y12, m14, c16, r18, g20, b22')
 if length(muv)>6
  figure, plot(r,AB(:,1:6)), hold on, plot(r,AB(:,7:length(muv)),'.-'),
 else
  figure, plot(r,AB),
 end
 title(' y0, m2, c4, r6, g8, b10, y12, m14, c16, r18, g20, b22')
 pausak
  figure, plot(muv(2:end),AB(1,2:end),'ro-'), 
  %a=axis;
  %ama=max(AB(1,2:end));
  %ami=min(AB(1,2:end));
  %a(3)=ami;
  %a(4)=ama;
  %axis(a)
  pausak
end

%figure, plot(muv,AB,'o-'), pausak
' sha_vortex fine', keyboard
