%ifp=10
%if istrumix==1
% aS=a;
% bS=b;
%
% a=aloc;
% b=bloc;
%end
%igint=0;
%' sha_scalini', keyboard

iplo=1;
iplo=0;

clear A B C D CD AB


an_veti= an_vetiVet(ndis);
an_vetu=an_vetuVet(ndis);
[an_veti' an_vetu']

  aa=aloc;
  ab=bloc;
  igint=1;
   ng=fix(301*Par.re/4);
  if igint==1

   [rv,wi]=gauleg(aa,ab,ng);
  else
   rv=linspace(aa,ab,ng);
  end
  r=rv/(kcav);  
  sgimp=ones(size(r));

muv=[0:2:2*nubesu];
muv=[0:pasnu:3*nubesu];
%muv=[0:2:18];
AB=ones(2,length(muv))*NaN;
CD=ones(2,length(muv))*NaN;
im=0;
ab0=an_vetu-an_veti;
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
 
 %'fine scalini', keyboard