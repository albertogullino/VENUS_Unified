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
    dmn=mu-nu+2;
    fim=find(dmn==muv);
    if length(fim)==1
     mfats=AB(:,fim);
     Bpm(:,jmu,jnu)=AB(:,fim);
     Bmp(:,jnu,jmu)=AB(:,fim);
     if ipolar==0
      D(:,jmu,jnu)=CD(:,fim);
     end
    else
%     disp('errore B mu in sha_oxi ')
%     pausak
    end
   end  %if   
%  end
 end
B=Bmp;

%'fineMPNB', keyboard