%' entr/pio modangle', keyboard

if ~exist('ipost')==1
 ipost=0;
end
clear besm besp
  npk=length(KK);
  vfa=ones(npk,1);
  nk1v=diff(ldap);

  if iLP==0
    iztm=1;
    mbv=[-1+nubesi:nubesu+1]+(mmsav-mm);
    mbv1=[-2+nubesi:nubesu+2]+(mmsav-mm);
    if mm~=mmsav
     if itetm==1
      pola=-1
     else
      pola=1
     end
       pimu=pimu-(mmsav-mm);
       KKt=[KKt; KKt];
       if ipost==0
        pes=[pes; pes];
        Pus=[Pus Pus+length(Pus)];
       end 
       Pusas=Pus;
       KKtz=[KKv*0; KKv];
       KKtz=KKtz(Pus);
       KKth=[KKv; KKv*0];
       KKth=KKth(Pus);
       lKAn=lKAn*2;

    end
    lbv=length(mbv);
%' entro modangle 1', keyboard
    nup1=numodi;
    nup1=(nubesu-nubesi)/pasnu+1;
    ibm=[1:pasnu:lbv-2];
    ibp=[3:pasnu:lbv];
    ibz=[2:pasnu:lbv-1];
%    if ipolar==0
     ibm1=[ibm -ibm];
     ibp1=[ibp -ibp];
     ibz1=[ibz -ibz];
%    else
%     ibm1=[ibm];
%     ibp1=[ibp];
%     ibz1=[ibz];
%    end

if ipolar==0
 irep=4;
else
 irep=2;
end 
   clear besp besm besz
    for ur=1:Nx
      x0=x(ur)*KK*kcav0;
      bes=besselJ(mbv,x0);
      besm0=reshape(bes(:,ibm),npk*nup1,1);
      besp0=reshape(bes(:,ibp),npk*nup1,1);
%      besp(:,ur)=[besp0; -besp0];
      besp(:,ur)=repmat(besp0,irep,1);
      besm(:,ur)=repmat(besm0,irep,1);
      if iztm==1
       besz0=reshape(bes(:,ibz),npk*nup1,1);
%       besz(:,ur)=[besz0*0; besz0];
       besz(:,ur)=repmat(besz0,irep,1);
       fkt=sqrt(KKtz.^2./(1-KKtz.^2));
       fkth=sqrt(KKth.^2./(1-KKth.^2));
      end
    end
   if iredmat==0
    Pusc=Pus;
   else
    Pusc=Pusas;
   end
    besp=besp(Pusc,:);
    besm=besm(Pusc,:);
    if iztm==1
     besz=besz(Pusc,:);
    end


    if pola==1
%     fmu=cos(mbv'*fian);
%     gmu=-sin(mbv'*fian);
     fmu=[cos(mbv'*fian); cos(mbv'*fian)];
     gmu=[-sin(mbv'*fian); -sin(mbv'*fian)];
    elseif pola==-1
%     fmu=sin(mbv'*fian);
%     gmu=cos(mbv'*fian);
     fmu=[sin(mbv'*fian); sin(mbv'*fian)];
     gmu=[cos(mbv'*fian); cos(mbv'*fian)];
    elseif pola==0
     fmu=[cos(mbv'*fian); sin(mbv'*fian)];
     gmu=[-sin(mbv'*fian); cos(mbv'*fian)];
    end
%     ' modangle ', keyboard

   Mvefm0=[];
   Mvegm0=[];
%   icoi=0;
   for nmoid=ibm1

    if nmoid<0
     nmoi=abs(nmoid)+length(mbv); 
     seg=segem;
    else 
     nmoi=abs(nmoid); 
     seg=1;
    end
    if ipolar==0
     seg=1;
    end
    Mdu=seg*vfa*fmu(nmoi,:);
    Mvefm0=[Mvefm0; Mdu];
    Mdu=seg*vfa*gmu(nmoi,:);
    Mvegm0=[Mvegm0; Mdu];
   end
%   Mvefm=diag(AnFs)*[ Mvefm0;  segem*Mvefm0];
%   Mvegm=diag(AnFs)*[ Mvegm0;  segem*Mvegm0];
%
%   Mvefp=diag(AnFs)*[ Mvefp0; -segem*Mvefp0];
%   Mvegp=diag(AnFs)*[-Mvegp0;  segem*Mvegp0];

   Mvefp0=[];
   Mvegp0=[];
%   icoi=0;
   for nmoid=ibp1
%    icoi=icoi+1;
%    vfa=ones(nk1v(icoi),1);
    if nmoid>0
     segf=1;
     segg=-1;
    else 
     segf=-segem;
     segg=segem;
    end
    if ipolar==0
     segf=1;
     segg=1;
    end
    if nmoid<0
     nmoi=abs(nmoid)+length(mbv); 
    else 
     nmoi=abs(nmoid); 
    end    
%    [nmoid nmoi], pausak
    Mdu=segf*vfa*fmu(nmoi,:);
    Mvefp0=[Mvefp0; Mdu];
    Mdu=segg*vfa*gmu(nmoi,:);
    Mvegp0=[Mvegp0; Mdu];
   end
%   Mvefm=diag(AnFs)*[ Mvefm0;  segem*Mvefm0];
%   Mvegm=diag(AnFs)*[ Mvegm0;  segem*Mvegm0];
%
%   Mvefp=diag(AnFs)*[ Mvefp0; -segem*Mvefp0];
%   Mvegp=diag(AnFs)*[-Mvegp0;  segem*Mvegp0];

   if iztm==1
     Mvez0=[];
     Mvhz0=[];
%     icoi=0;
     for nmoid=ibz1
         if nmoid<0
          nmoi=abs(nmoid)+length(mbv); 
          seg=segem;
         else 
          nmoi=abs(nmoid); 
          seg=1;
         end
%      icoi=icoi+1;
%      vfa=ones(nk1v(icoi),1);
      Mdu=segem*vfa*fmu(nmoi,:);
      Mvez0=[Mvez0; Mdu];
      Mdu=vfa*gmu(nmoi,:);
      Mvhz0=[Mvhz0; Mdu];
     end
   end
%'modal', keyboard
   lPu=length(Pusc)/2;
%   Pus1=Pusc(1:lPu);
%   if ipolar==0 & iLP==0
    Pus1=Pusc;
    Pus2=Pus1;
    
%   else
%    Pus2=Pus1+Pus1(end); 
%   end
  
%   lpp=find(Pusc==120);
%   Pus1=Pusc(1:lpp);
%   Pus2=Pusc(lpp+1:end)-120;
%' modan ', keyboard
%' modan ', keyboard
   if ipolar==0 & iLP==0
    Mvefm0s= Mvefm0; 
    Mvegm0s=Mvegm0; 
    Mvefp0s=Mvefp0;    
    Mvegp0s=Mvegp0;   

%   Mvefm=diag(AnFs)*[ Mvefm0;  segem*Mvefm0];
%   Mvegm=diag(AnFs)*[ Mvegm0;  segem*Mvegm0];
%
%   Mvefp=diag(AnFs)*[ Mvefp0; -segem*Mvefp0];
%   Mvegp=diag(AnFs)*[-Mvegp0;  segem*Mvegp0];
    pu1=1:size(Mvefp0,1)/2;
    pu2=size(Mvefp0,1)/2+1:size(Mvefp0,1);

    Mvefm0=[Mvefm0(pu1,:); segem*Mvefm0(pu1,:); Mvefm0(pu2,:); segem*Mvefm0(pu2,:)]; 
    Mvegm0=[Mvegm0(pu1,:); segem*Mvegm0(pu1,:); Mvegm0(pu2,:); segem*Mvegm0(pu2,:)]; 
    Mvefp0=[Mvefp0(pu1,:); -segem*Mvefp0(pu1,:); Mvefp0(pu2,:); -segem*Mvefp0(pu2,:)]; 
    Mvegp0=[-Mvegp0(pu1,:); segem*Mvegp0(pu1,:); -Mvegp0(pu2,:); segem*Mvegp0(pu2,:)]; 


    if iztm==1
     Mvez0s=Mvez0;   
     Mvhz0s=Mvhz0;   
     Mvez0=[Mvez0(pu1,:); Mvez0(pu1,:); Mvez0(pu2,:); Mvez0(pu2,:)];
     Mvhz0=[Mvhz0(pu1,:); Mvhz0(pu1,:); Mvhz0(pu2,:); Mvhz0(pu2,:)];
    end  

   end

    Mvefm0su= Mvefm0; 
    Mvegm0su=Mvegm0; 
    Mvefp0su=Mvefp0;    
    Mvegp0su=Mvegp0;   

   Mvefp0=Mvefp0(Pus1,:);
   Mvegp0=Mvegp0(Pus1,:);

   Mvefm0=Mvefm0(Pus2,:);
   Mvegm0=Mvegm0(Pus2,:);

   
   if iztm==1
    Mvez0=Mvez0(Pus1,:);
    Mvhz0=Mvhz0(Pus1,:);
   end

%   'Mvef', keyboard
  else  %iLP=1

    ibm=[1:pasnu:lbv];
    nup1=length(ibm);;
    mbvc=mbv(ibm);

    for ur=1:Nx
      x0=x(ur)*KK*kcav0;
      bes=besselJ(mbvc,x0);
      besm0=reshape(bes,npk*nup1,1);
      besm(:,ur)=[besm0];
    end


   if pola==1
    fmu=cos(mbvc'*fian);
   elseif pola==-1
    fmu=sin(mbvc'*fian);
   elseif pola==0
    fmulo=[cos(mbvc'*fian); sin(mbvc'*fian)];
    fmu=[cos(mbvc'*fian); sin(mbvc(2:end)'*fian)];
   end
   if pola~=0
    fmulo=fmu;
   end
%' fmu', keyboard
   Mvefm0=[];
%   icoi=0;
%   for nmoi=1:length(ibm)
   for nmoi=1:min(size(fmulo))
%    icoi=icoi+1;
%    vfa=ones(nk1v(icoi),1);
    Mdu=vfa*fmulo(nmoi,:);
    Mvefm0=[Mvefm0; Mdu];
   end

   Pus0=Pus;
%   if ipolar==0
%    Pus=[Pus0 Pus0+length(Pus0)];
%   else
   if iredmat==0
    Pusc=Pus0;
   else
    Pusc=Pusas;
   end
%   end

   Mvefm0=Mvefm0(Pusc,:);
   if ipolar==0
    besm=[besm; besm];
   end
   besm=besm(Pusc,:);

  end  %iLP

  if isi==1
   Fs=sqrt(KKt(Pus));
  else
   Fs=1;
  end

%'angle', keyboard