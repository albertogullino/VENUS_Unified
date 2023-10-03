if ~exist('iallo')
 iallo=2; 
 if iLP==1
  iallo=1;
 end
end 

iploloc=0;
lmbv=length(mbv);
%' po_new passo', keyboard
nk1max=length(KK);
if iredmat==1
 Pusff=Pusas;
else
 Pusff=Pus;
end
te0g=90;
%Aep=Afz(:,iFFsect)+Afzf(:,iFFsect);
%Ahp=(-Afz(:,iFFsect)+Afzf(:,iFFsect));
Aep=Andu;
Ahp=AnduH;
   if ired_ret==1
    A1=zeros(ldap(end),1);
    A1(Pured1)=Aep;
    A2=zeros(ldap(end),1);
    A2(Pured1)=Ahp;
    Aep=A1;
    Ahp=A2;
   end

    if mm~=mmsav0
%    if mm~=mmsav
     if itetm==1
      Aep=[Aep; zeros(size(Aep))];
      Ahp=[Ahp; zeros(size(Ahp))];
     else
      Aep=[zeros(size(Aep)); Aep];
      Ahp=[zeros(size(Ahp)); Ahp];
     end
    end


%Ahp=segdir*(Afz(:,iFFsect)-Afzf(:,iFFsect));


if iant==0
 Anf=Aep.*pes(Pus);
 Anh=Ahp.*pes(Pus);
else
 Anf=Aep;
 Anh=Ahp;
end


fkttn=KKt./sqrt(1-(rrfa*KKt).^2);



lim=1-(rrfa*KKv).^2;
kacce=find(lim<0);
kaccea=[kacce; kacce+length(lim)];

 if ipolar==0
  irep1=2;
 else 
  irep1=1;
 end
if iLP==0
SQ=sqrt(1-(KKv).^2);

Znorn=[1./SQ; SQ];
Znorn=repmat(Znorn,irep1,1);
Ynorn=rr./Znorn;
%Ynorn(kaccea)=0;
Ynorn=Ynorn(Pusff);

SQi=repmat(SQ,irep1*2,1);
fkttn=KKt./SQi;

%fkttn(kaccea)=0;
fkttn=fkttn(Pusff);

 Anfz=-segem*2*(-j)*Anh.*fkttn*fatEz;
 AnfHz=2*(-j)*Anf.*fkttn;
 AnfH=Anh.*Ynorn;
else
Ynorn=rr;
fkttn=KKt;
fkttn=repmat(fkttn,irep1,1);
end


%fkttn=fkttn(Pusff);



if iFF==1
%      z0c=input(' distanza ff in cm = ');
      z0=z0c*1e4;
      lambdaff=lambda*(1-ze)/rff;
      k0=2*pi/lambdaff;
      teRg=asin(Kff*rrfa)*180/pi;
      fik=find(teRg>te0g);
      if length(fik)>0
       fit=fik;
       Anf(fit)=0;
       figure, plot(abs(Anout)), pausak
      end
     if iFFte==0
% vecchio ff
      rox=z0c*tan(teR);
      R=z0./cos(teR);
%      Fdu=-1./(R.*tan(teR)).*exp(-j*k0*R+j*pi/4);
%      Fd0=cos(teR)./(2*pi*lambdaff*R).*exp(-j*k0*R+j*pi/4);
      Fd0=-cos(teR)./(4*lambdaff*R);  % the exponential factor cancel out in Point
      Fd=diag(Fd0);
     else
% nuovo ff
      Rd=z0;
%      Fd0=1/(2*pi*lambdaff*Rd)*exp(-j*k0*Rd+j*pi/4);
      Fd0=1./(4*lambdaff*Rd);  % the exponential factor cancel out in Point
      Fd=Fd0*diag(cos(teR));
     end
      Fdz=diag(tan(teR));

      lK=length(KK);
      lKv=diff(ldap);
      ldap1=ldap;
      ldap1=ldap;
%      ldap1(end)=ldap1(end)-1;
      ldap1=ldap1+1;
      Pus1=Pus;
      Pus1=[Pus1 Pus1(end)+1];

    if iLP==0

      Efx=0;
      Efy=0;
      Efz=0;
      Hfx=0;
      Hfy=0;
      Hfz=0;


      AnfHz=2*(-j)*Anf.*fkttn;


      lkr=1;
      if exist('krv')
        lkr=length(krv);
      end
      lkr1=1;
      Puf=1:nk1max;
      Pudu=repmat(Puf,1,numodi*iallo*iall)';
%      for imi=1:length(ibm)*lkr



%      for imi=1:(lbv-1)/pasnu*lkr
%' attenzione modificato con numodi !', 
%pausak
ibp2=abs(ibp1);
ibm2=abs(ibm1);
ibz2=abs(ibz1);
%      for imi=1:length(mbvero)
motut=fix(abs(Moe(ldap(2:end))));
mmlocvec=mbv(ibz2);       

   for imi=1:length(motut);
       imiver=ceil(imi/lkr);
       imi1=imiver-1;
       mmloc=motut(imiver);
       mmfi=find(mmloc==mmlocvec);
       imivero=mmfi(1);
       mmp=ibp2(imivero);
       mmm=ibm2(imivero);
       mmi=ibp2(imivero)-1;
       nui=mmloc;
       num=mmloc-1;
       nup=mmloc+1;
%       lK=lKv(imi);
%       pumi=ldap(imi)+1:ldap(imi+1);
       lK=lKv(imiver);
       pumi=ldap(imiver)+1:ldap(imiver+1);
       pumif=Pus1(ldap1(imi)):Pus1(ldap1(imi+1)-1);
%'cont ', keyboard
      nuv_loca=fix(abs(Moe(pumi(1))));
      nuv_locr=(real(Moe(pumi(1))));
      nuv_loci=(imag(Moe(pumi(1))));

%     '[imi, mmloc], keyboard '
%     [imi, mmloc], 
%     keyboard
     if  mmloc==nuv_loca
       AnEd=Anf(pumi);
       AnE=zeros(nk1max,1);
       AnE(Pudu(Pusff(pumi)))=AnEd;
       
%       AnMd=Anf(pumi+lKAn/2);
%       AnM=zeros(nk1max,1);
%       AnM(Pudu(Pusff(pumi)))=AnMd;



       Imm=Fd*reshape(Im(:,:,mmm),npkf,npkF);
       Imp=Fd*reshape(Im(:,:,mmp),npkf,npkF);
       Imi=Fd*reshape(Im(:,:,mmi),npkf,npkF);              
       if  nuv_locr<0 | nuv_loci<0
        sed=-segem;
        ses=segem;
        stetm=0;
       elseif nuv_locr>0 | nuv_loci>0
        sed=1;
        ses=1;
        stetm=1;
       end
       
       if abs(nuv_loci)>0
        stfmu=lmbv;
       else
        stfmu=0;
       end
       mmp=mmp+stfmu;
       mmm=mmm+stfmu;
       mmi=mmi+stfmu;
       
%  Campo elettrico       
%       AD=(AnE-segem*AnM);
%       AS=(AnE+segem*AnM);
       AD=sed*AnE;
       AS=ses*AnE;
       Efx=Efx+((Imp*AD.*j^nup)*fmu(mmp,:)+(Imm*AS.*j^num)*fmu(mmm,:));
       Efy=Efy+((-Imp*AD.*j^nup)*gmu(mmp,:)+(Imm*AS.*j^num)*gmu(mmm,:));

       if stetm==0
%       AnMd=Anfz(pumi+lKAn/2);
        AnMd=Anfz(pumi);
        AZ=zeros(nk1max,1);
        AZ(Pudu(Pusff(pumi)))=AnMd;
        if iploloc==1
         figure, plot(abs(AnE),'r.'),  hold on,  %plot(abs(AnM),'g.'),
         plot(abs(AZ),'c.'),
         title(' Campo elettrico E M Z'), pausak
        end

        Efz=Efz+(Imi*AZ)*j^nui*fmu(mmi,:);
       end  % stetm

%  Campo magnetico

       AnEd=AnfH(pumi);
%       AnMd=AnfH(pumi+lKAn/2);
       AnE=zeros(nk1max,1);
       AnE(Pudu(Pusff(pumi)))=AnEd;
%       AnM=zeros(nk1max,1);
%       AnM(Pudu(Pusff(pumi)))=AnMd;

       AD=sed*AnE;
       AS=ses*AnE;
       Hfy=Hfy+((Imp*AD.*j^nup)*fmu(mmp,:)+(Imm*AS.*j^num)*fmu(mmm,:));
       Hfx=Hfx-((-Imp*AD.*j^nup)*gmu(mmp,:)+(Imm*AS.*j^num)*gmu(mmm,:));

       if stetm==1 
        AnMd=AnfHz(pumi);
        AZ=zeros(nk1max,1);
        AZ(Pudu(Pusff(pumi)))=AnMd;

        if iploloc==1
         figure, plot(abs(AnE),'r.'), hold on,  %plot(abs(AnM),'g.'),
         plot(abs(AZ),'c.'),
         title(' Campo magnetico '), pausak
        end
        Hfz=Hfz+(Imi*AZ)*j^nui*gmu(mmi,:);
       end % stetm  
      end

     end



      Ef=(Efx).^2+(Efy).^2+(Efz).^2;
%      Ef=(Efx).^2+(Efy).^2;
%      Ef=(Efy).^2;
%      'passo', keyboard
      mfxy=(max(max(abs(Ef))));
      mfxyc=sqrt(max(max((Ef))));
      Ef=sqrt(abs(Ef/mfxy));

    else   %iLP



      Efx=0;
      Hfy=0;

      ibm=[1:pasnu:lbv];
      mbvc=mbv(ibm);
      if ipolar==0
       ibm=[ibm ibm(2:end)];
       mbvc=[mbvc mbvc(2:end)];
      end


      lPus=0;
%      for imi=1:length(ibm)
      lkr=1;
      if exist('krv')
        lkr=length(krv);
      end
      lkr1=1;
      Puf=1:nk1max;
      Pudu=repmat(Puf,1,numodi)';
      for imi=1:length(ibm)*lkr
%      for imi=1:lkr:length(ibm)*lkr
%      for imi=1:lkr:length(ibm)*lkr
       imiver=ceil(imi/lkr);
       mmloc=mbvc(imiver);

       mmi=ibm(imiver);
       nui=mmloc;
       lK=lKv(imi);
       pumi=ldap(imi)+1:ldap(imi+lkr1);
       pumif=Pus1(ldap1(imi)):Pus1(ldap1(imi+lkr1)-1);

       AnEd=Anf(pumi);
%       figure, plot(abs(Anf)), hold on, plot(pumi,abs(AnEd),'r.'),  pausak
%       AnE=[zeros(nk1max-length(AnEd),1); AnEd];
       AnE=zeros(nk1max,1);
       AnE(Pudu(Pusff(pumi)))=AnEd;
%' imi  mmi nui', 
%[ imi  mmi nui], 
%       figure, plot(abs(AnE),'r'); pausak

       Imi=Fd*reshape(Im(:,:,mmi),npkf,npkF);

       AD=AnE;

%'Efx', keyboard
       Efx=Efx+(Imi*AD.*j^nui)*fmu(imiver,:);
%       Efx=Efx+(Imi*AD)*fmu(imiver,:);
%       Efxm{imiver}=(Imi*AD.*j^nui)*fmu(imiver,:);

%       figure, plot(abs(Efx)), pausak


%       AnEd=AnfH(pumi);
%       AnE=[zeros(nk1max-length(AnEd),1); AnEd];
       AnEd=Anh(pumi);
       AnE=zeros(nk1max,1);
       AnE(Pudu(Pusff(pumi)))=AnEd;

       AD=AnE;

       Hfy=Hfy+(Imi*AD.*j^nui)*fmu(imiver,:);
%       Hfy=Hfy+(Imi*AD)*fmu(imiver,:);
%       ' po_pro', pausak

      end


    end

%'z0'
%keyboard
isal=1;
if isal==0
  fx=find(imag(roxt)==0);
  Xdu=roxt(fx)*fm0*1e-4;
  Ydu=roxt(fx)*gm0*1e-4;
  Ef=Ef(fx,:);

  if iFFte==1
   X=atan(Xdu/z0c)*180/pi;
   Y=atan(Ydu/z0c)*180/pi;
   axli=60;
  else
   X=Xdu;
   Y=Ydu;
   axli=z0c;
  end

   fiazi=atan(sqrt(Xdu.^2+Ydu.^2)/z0c);
   fipla=ones(size(fx))*fian0;
else

  if iFFte==1

   axli=60;
  else
   axli=z0c;
  end

end %isal

% Pvectz=(Ex.*conj(Hy)-Ey.*conj(Hx));
% Pvectx=(Ey.*conj(Hz)-Ez.*conj(Hy));
% Pvecty=-((Ex.*conj(Hz)-Ez.*conj(Hx)));


     if iLP==0
      Pvectx=real(Efy.*conj(Hfz)-Efz.*conj(Hfy));
      Pvecty=-(real(Efx.*conj(Hfz)-Efz.*conj(Hfx)));
      Pvectz=real(Efx.*conj(Hfy)-Efy.*conj(Hfx));
      Pointing=sqrt(Pvectx.^2+Pvecty.^2+Pvectz.^2);

     else
%      Pvectz=(real(Efx.*conj(Hfy)));
%      '    Attenzione: cambio in po_pro '
      Pvectz=abs(Efx.*conj(Hfy));
     end
%     Ppo=Pointing/Pointing(1,1);
%     Ppo=Pvectz/Pvectz(1,1);

% figure, plot(X1(:,1),abs(Efx(:,46)/Efx(1,1)).^2,X1(:,1),...
% abs(Ppo(:,46)))
%load sa1
% hold on, plot(X1(:,1),Ec,'.',X1(:,1),Pc,'.')

% Ec=abs(Efx(:,46)/Efx(1,1)).^2;
% Pc=abs(Pvectz(:,46)/Pvectz(1,1));
% save sa Pc Ec

iskf=1;
if iskf==0
  figure
  subplot(1,3,1)
  surf(X,Y,(abs(Efx)/mfxy).^2),
  shading('interp'), view(0,90),
  axis square, axis equal, grid,
  axis([-1 1 -1 1]*axli/2),

  subplot(1,3,2)
  surf(X,Y,(abs(Efy)/mfxy).^2),
  shading('interp'), view(0,90),
  axis square, axis equal, grid,
  axis([-1 1 -1 1]*axli/2),

%  figure
%  surf(Xte,Yte,(abs(Efy)/mfxy).^2),
%  shading('interp'), view(0,90),
%  axis square, axis equal, grid,

  subplot(1,3,3)
  surf(X,Y,(abs(Ef)/mfxy).^2),
  shading('interp'), view(0,90),
  axis square, axis equal, grid,
  axis([-1 1 -1 1]*axli/2),
  xlabel(' x cm '), ylabel(' y  cm'),
  drawnow
end   %iskf
 if ifp>=-3
  pausak
 end
end %iFF

%'end po_new' , keyboard
