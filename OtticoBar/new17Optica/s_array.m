
igintsa=igint;
igint=0;
ifalso=-1;
icirc0=0;
clear A B CD AB

sgimp=1;
if istrumix==1
 avero=aloc/kcav0;
end


npfia=11*numodiacc;
npfia=51;
%'!!!!!!!!!!!!!!!!! in s_arraya:  npfia=11'
%npfia=11;
np00=31;

fi=linspace(0,2*pi,4*npfia+1)';
fi00=linspace(0,2*pi,8*npfia+1)';

%'QUI'
%keyboard

sh_type=Pshi{1};
Rx=Pshi{2}(1);
Ry=Pshi{3}(1);
if length(find(abs(sh_type)==1))==length(sh_type)
 sarray_c
% sarray_g
else
 sarray_g
end



muv=[0:2:2*(nubesu)];
%muv=[0:2:30];
AB=zeros(length(r),length(muv));
CD=AB;
im=0;
for mu=muv
im=im+1;
 for k=1:length(r)
  fid=fasid(k,:);
  fud=fasud(k,:);
  fini=find(isnan(fid)==0);
  finu=find(isnan(fud)==0);
  if length(fini)~=0 & length(fini)==length(finu)
   if mu==0
    Ad=sum(fud(finu)-fid(fini));
    Cd=0;
   else
    Ad=sum(sin(mu*fud(finu))-sin(mu*fid(fini)))/mu;
    Cd=-sum(cos(mu*fud(finu))-cos(mu*fid(fini)))/mu;
   end
%   if abs(Cd)>1e-6
%    Cd
%    pausak
%   end
   AB(k,im)=Ad/2;
   CD(k,im)=Cd/2;
  end
 end
 mfiA=mean(abs(AB(:,im)));
 if mfiA<2e-3
  AB(:,im)=0;
  ABac(im)=0;
 else
  ABac(im)=1;
 end
% fiA=find(abs(CD(:,im))<1e-14);
% CD(fiA,im)=CD(fiA,im)*0;
 mfiA=mean(abs(CD(:,im)));
 if mfiA<1e-4
  CD(:,im)=0;
 end

% fiA=find(abs(AB(:,im))<1e-7);
% if length(fiA)>ng/2
%  AB(:,im)=AB(:,im)*0;
% end
% fiA=find(abs(CD(:,im))<1e-7);
% if length(fiA)>ng/2
%  CD(:,im)=CD(:,im)*0;
% end
end
if ifp>1
 if length(muv)>6
  figure, plot(r,AB(:,1:6)), hold on, plot(r,AB(:,7:length(muv)),'.-'),
 else
  figure, plot(r,AB),
 end
 title(' AB: y0, m2, c4, r6, g8, b10, y12, m14, c16, r18, g20, b22')
 pausak
 if length(muv)>6
  figure, plot(r,CD(:,1:6)), hold on, plot(r,CD(:,7:length(muv)),'.-'),
 else
  figure, plot(r,CD),
 end
 title(' CD: y0, m2, c4, r6, g8, b10, y12, m14, c16, r18, g20, b22')
 pausak
end

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
%figure, plot(ro,AB,ro,CD)

igint=igintsa;
clear Mlo Flo Rsa Rlo
