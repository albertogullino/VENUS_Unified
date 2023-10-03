%ifp=10
%if istrumix==1
% aS=a;
% bS=b;
%
% a=aloc;
% b=bloc;
%end
%igint=0;
' sha_vortex', keyboard
iplo=0;
%Par=Ps.Par;
%mazim=Par.mazim;
%SogDen=Par.SogDen;

%mazim=4;
%SogDen=.2;

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

if ifp>0
firad=linspace(0,2*pi,501);
figure, plot(firad,cos(mazim*firad),an_veti,cos(mazim*an_veti),'go',an_vetu,cos(mazim*an_vetu),'ro'), 
hold on
plot(firad,ones(size(firad))*SogDen,'w--',firad,-ones(size(firad))*SogDen,'w--')
pausak
figure, plot(firad,sin(mazim*firad),an_veti,sin(mazim*an_veti),'go',an_vetu,sin(mazim*an_vetu),'ro'), 
hold on
plot(firad,ones(size(firad))*SogDen,'w--',firad,-ones(size(firad))*SogDen,'w--')
pausak
end
rv=[adis bdis];

muv=[0:2:2*nubesu];
%muv=[0:2:18];
AB=ones(2,length(muv))*NaN;
im=0;
ab0=mazim*(an_vetu(2)-an_veti(2));
for mu=muv
im=im+1;
 if mu==0
  AB(1:2,im)=ab0;
 else
  AB(1:2,im)=1/(2*mu)*sum(sin(mu*an_vetu)-sin(mu*an_veti));
 end
 %(sin(mu*an_vetu)-sin(mu*an_veti))', pausak
end

fi=find(abs(AB)<1e-10);
AB(fi)=0;

'qui',  keyboard

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

figure, plot(muv,AB,'o'), pausak
' sha_vortex fine', keyboard
