%if ifp==-10
% 'entro fie_newARM', keyboard
%end

if exist('PV')==1
 if PV{1}.D<1.5
  Cug.x=[];
  Cug.y=[];
  Cug.z=[];
 end
end 
iLP=iLP1;
%iLP=iLPr;
iint=1;
axlim=alim*180/pi*rr;  %1.5 fattore arbitrario
if axlim>40
 axlim=40;
end
ibaru=1;
rtetm=0;
nuazi=0;
mrad=0;
polca=0;
polratio=0;
saou=0;
iEz=1;
if iLP==1
 iEz=0;
end
fatEz=(rr/rfu)^2;
%'fienew', keyboard
imod=imod+1;
if ifp~=-4
 imod
end
iaoff=0;
iskim=0;        % serve per LP e modo fondamentale tipo sin(0*fi)
icut=0;
if ~exist('icaplo')
 icaplo=2;
end
if length(icaplo)==0
 icaplo=2;
end
if ~exist('Cug')
 Cug.x=0;
 Cug.y=0;
 Cug.z=0;
end
%  vgconv=3e10/rr;
  glosout=2*gg0;
  Gsimmi=gg0;
  fris=ze;
%  freq=fris;
  ADoutis=An0;
  ADouti=An0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% decido quali punti prendere dell'autovettore (puAc)

    taglio

% FINE decido quali punti prendere dell'autovettore
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calcolo tutti i campi nelle sezioni che servono

%' calcnafo', keyboard
if ~exist('gplam')
 gplam=0;
end

 calc_new_arm

%  FINE calcolo tutti i campi nelle sezioni che servono
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% assegno i vari parametri modali (se il modo e` buono: iskim=0)
%' sono qui ub fienew', keyboard

if iskim==0

   glostotal=glosout;


 if length(if_only_out)==0
   ipu_sav=2;
   ADc(1:length(Afz(:,2)),imod)=Afz(:,2)+Afzf(:,2);
   ADout(1:length(Afz(:,2)),imod)=Afz(:,iFFsect)+Afzf(:,iFFsect);
 else
   ipu_sav=1;
   icaplo=1;
   ADc(1:length(Afz),imod)=Afz+Afzf;
   ADout(1:length(Afz),imod)=Afz+Afzf;
 end
   %Gsov(imod)=glostotal+gplam*vgconv;
   Gsov(imod)=glostotal;
   Fsov(imod)=ze;
   IndRad(imod)=mm;
   Gsosa=Gsov;
   Fsosa=Fsov;
   
   
%   'qui', keyboard
iM2=1;
M2=1;
if iM2==1
 xri=xvero;
 xdx2=xri.^3*diff(xri(1:2));
 xdx=xri.*diff(xri(1:2));
 Eca=mean(abs(Etx(:,:,2).^2)+abs(Ety(:,:,2).^2),2);
 Es=xdx2*Eca;
 En=xdx*Eca;
 sig2=Es./En;

 z=10e3;
 tet=X(:,1)/180*pi;
 xf=z*tan(tet)';
 dx=diff(xf);
 dx=[dx(1) dx];
 xdx2f=xf.^3.*dx;
 xdxf=xf.*dx;
 Ecaf=mean(abs(Ef).^2,2);
 Esf=xdx2f*Ecaf;
 Enf=xdxf*Ecaf;
 sig2f=Esf./Enf;
 
 tet=X(:,1)/180*pi;
 dtet=tet.*[0; diff(tet)];

 dsi0=2*sqrt(2*sig2); 
 dsi0f=2*sqrt(2*sig2f);

 M2=pi/(4*lambda*z)*dsi0.*dsi0f;
 if ifp==-10
%   ' qui per M2 ', keyboard
%   ' qui per M2 ', keyboard
 end  
  M2v(imod)=M2;

end  

   Nazim(imod)=nuazi;
   Nrad(imod)=mrad;

%   Nazim(imod)=1;
%   Nrad(imod)=1;
%   ' Nrad ', keyboard

   tyE(imod)=1;

   if idyn>=1
    TL_ef
    W=0;
   end

   if i2D==3

     s=size(Etot);
     p1=1:s(1);
     p2=1:s(2);
     EDc.x(p1,p2,imod)=Etx(:,:,ipu_sav);
     EDout.x(p1,p2,imod)=Etx(:,:,1);
     E2xo=Etx(:,:,iFFsect);
     E2xp=Etx(:,:,icaplo);
     E2zo=Etz(:,:,iFFsect);


     if iLP==0
      EDc.y(p1,p2,imod)=Ety(:,:,ipu_sav);
      EDout.y(p1,p2,imod)=Ety(:,:,1);
      E2yo=Ety(:,:,iFFsect);
      E2zo=Etz(:,:,iFFsect);
      E2zp=Etz(:,:,icaplo);
      E2yp=Ety(:,:,icaplo);
     end

     if iLP==0
      if max(max(real(E2y)))>max(max(real(E2x)))
       tyE(imod)=-1;
      end
     end

   else

     EDc.x(:,imod)=Etx(:,ipu_sav);
     EDout.x(:,imod)=Etx(:,1);
     if iLP==0
      EDc.y(:,imod)=Ety(:,ipu_sav);
      EDout.y(:,imod)=Ety(:,1);
     end

   end
%'fienoiw vedp', keyboard
      Ap=abs(Anout)/max(abs(Anout));
      ApQ=abs(AnQW)/max(abs(AnQW));

%    if ifp<=-10 | ifp==-4
iskv=0;

   if iskv==0
    aax=0.8*max(xvero);

%' propio qu', keyboard

    if ifp<=-10
    ifps=ifp;
    ifp=1;
     fgsav=figure;

      if iLP==0
       nsub=3
       pos=20;
      else
       nsub=2
       pos=20;
      end
      if pola==1
       set(gcf,'Position',pograp);
       pograp(1)=pograp(1)+pos*sinc;
       if sinc<1
        if pograp(1)<50
         pograp(1)=950;
         pograp(2)=10;
        end
       else
       if pograp(1)>950
        pograp(1)=50;
        pograp(2)=10;
       end
       end
      else
       set(gcf,'Position',pogram);
       pogram(1)=pogram(1)+pos*sinc;
      end

      subplot(nsub,nsub2,1)



      if pola==1
       plot(Ap.*abs(Ofin),'r.-'), hold on, plot(ApQ,'w') ;
      else
       plot(Ap.*abs(Ofin),'g.-'), hold on, plot(ApQ,'w') ;
      end
%   figure,       plot(Ap,'g.-'), hold on, plot(ApQ,'w'),        plot(puaud,Ap(puini),'co'),
      a=axis;
      if (a(2)<10)
       a(2)=100;
       axis(a);
       AX=a;
      end
      if ~exist('kretvero')
       kretvero=0;
      end
%if kretvero==0
%      plot(KKvd(Pusas),'.');
%else
%     if ired_ret==1
%      plot(KKt(Pusas(Pured1)),'.');
%     else 
%      plot(KKt(Pusas),'.');
%     end 
%end
       puau=1:length(Ap);
       suee=sum(sum(Ef));
       
%   figure,       plot(Ap,'g.-'), hold on, plot(ApQ,'w'),        plot(puaud,Ap(puini),'co'),
       if suee>0 & ired_ret~=1
       puini=ldap1(1:end-1);
       puaud=puau(puini);
       plot(puaud,Ap(puini),'co'),
       end
      if iLP==0
       shi=0;
       shi1=-.08;
      else
       shi=-.2;
       shi1=0;
      end
%      stri=[' gain = ',num2str(glosout,'%0.4e')];
%      text(fix(length(Anout)/2.5),1.27+shi,stri);
%      stri=[' wav = ',num2str(lambda/(1+ze),'%0.4e')];
%      text(fix(length(Anout)/2.5),1.17+shi,stri);
%      stri=[' freq = ',num2str(ze,'%0.4e')];
%      text(fix(length(Anout)/2.5),1.07+shi,stri);
      zel=ze*lambda*1000;
      if ifiez<2
       glosp=glosout/vgconv;
      else
       glosp=2*ggg1/vgconv;
      end
      glosp=glosp;
      stri=[' gain = ',num2str(glosp,'%0.4e'),',   Dlam = ',num2str(zel,'%0.4e')];
      text(-.25,1.18-shi1,stri);
%      stri=[' alim = ',num2str(alim,'%0.2f'),',  Nk = ',num2str(nK_dis,'%3i'),', N-acc = ',num2str(numodiacc+1,'%2i')];
%      stri=[stri,' gain per QW = ',num2str(glosp/(NQW*fatqw),'%0.4e')];
      laver=lambda*1000/(1+ze);
      stri=['lam= ',num2str(laver,'%0.4e'),',     alim = ',num2str(alim,'%0.2f'),',  Nk = ',num2str(nK_dis,'%3i'),', N-acc = ',num2str(numodiacc+1,'%2i')];
%      stri=[stri,' gain per QW = ',num2str(glosp/(NQW*fatqw),'%0.4e')];
      text(-.25,1.08,stri);
      %'tedt', keyboard
      ax=axis;
      ax(2)=length(Anout+5);
      axis(ax)
      if exist('AX')
       axis(AX)
      end
%      text(fix(length(Anout)/2.5),1.27,stri);
%      text(fix(length(Anout)/2.5),1.27+shi,stri);
%      stri=[' wav = ',num2str(lambda/(1+ze),'%0.4e')];
%      text(fix(length(Anout)/2.5),1.17+shi,stri);
%      stri=[' freq = ',num2str(ze,'%0.4e')];
%      text(fix(length(Anout)/2.5),1.07+shi,stri);


       titl=' E_x  Output';
       nplo=2;
       subplot(nsub,nsub2,nplo+(nsub2-1)*(nplo-1)+1)
%       subplot(nsub,nsub2,3)
       ibar=0;

       if i2D==3
        map_fnew(XP,YP,E2xo,aax,Cug.x,Cug.y,Cug.z,titl,ibaru,iaoff)


% figure, map_fnew(XP,YP,Ef,aax,Cug.x,Cug.y,Cug.z,'forward',ibar,iaoff)
% figure, map_fnew(XP,YP,Eb,aax,Cug.x,Cug.y,Cug.z,'backward',ibar,iaoff)


%        figure
%        map_fnew(XP,YP,Exdu,aax,Cug.x,Cug.y,Cug.z,titl,ibar,iaoff)

%        figure
%        map_fnew(XP,YP,E2xo.^2,aax,Cug.x,Cug.y,Cug.z,titl,ibar,iaoff)

%        figure
%        map_fnew(XP,YP,E2yp,aax,Cug.x,Cug.y,Cug.z,titl,ibar,iaoff)

%        figure
%        map_fnew(XP,YP,real(Exdu),aax,Cug.x,Cug.y,Cug.z,titl,ibar,iaoff)
       else
        plot(xvero,Etcf)
        title(titl)
       end

       if iEz==1
        mXe=max(max(E2xp));
        mYe=max(max(E2yp));
        if mXe>mYe
         E2mp=E2xp;
         titlm=' E_x  QW';
        else
         E2mp=E2yp;
         titlm=' E_y  QW';
        end
       end

       if iEz==0
        if nsub2==2
         nplo=2;
         subplot(nsub,nsub2,nplo+(nsub2-1)*(nplo-1))
         ibar=0;
         iaoff=0;
         if i2D==3
          if iint==2
           Eint=abs(E2xp).^iint;
           titl=' |E_x|^2  QW';
          else
           Eint=E2xp;
           titl=' E_x  QW';
          end
          map_fnew(XP,YP,Eint,aax,Cug.x,Cug.y,Cug.z,titl,ibar,iaoff)
%          iaoff=1;
         else
          plot(xvero,Etcf)
          title(titl)
         end
        end
%        keyboard
       else
        if nsub2==2
         titl=titlm;
         nplo=2;
         subplot(nsub,nsub2,nplo+(nsub2-1)*(nplo-1))
         ibar=0;
         iaoff=0;
         if i2D==3
          if iint==2
           Eint=abs(E2mp).^iint;
          else
           Eint=E2mp;
          end
          map_fnew(XP,YP,Eint,aax,Cug.x,Cug.y,Cug.z,titl,ibar,iaoff)
%          iaoff=1;
         end
        end

       end

       subplot(nsub,nsub2,2)
       if iEz==1
%        titl='                       E_z  Output';
        titl='';
        aax=0.8*max(xvero);
        map_fnew(XP,YP,E2zo,aax,Cug.x,Cug.y,Cug.z,titl,ibaru,iaoff)
%        'qui Ez', keyboard
       end


       if iLP==0 & i2D==3
%        titl=' |E|^2 ';
%        subplot(nsub,1,3)
        titl=' E_y Output';
        nplo=3;
        subplot(nsub,nsub2,nplo+(nsub2-1)*(nplo-1)+1)
%        subplot(nsub,nsub2,5)
         map_fnew(XP,YP,E2yo,aax,Cug.x,Cug.y,Cug.z,titl,ibaru,iaoff)


        if iEz==0
        subplot(nsub,nsub2,nplo+(nsub2-1)*(nplo-1))
         if nsub2==2
          titl=' E_y  QW';
          nplo=3;
          ibar=0;
          if i2D==3
           aax=0.8*max(xvero);
           map_fnew(XP,YP,E2yp.^iint,aax,Cug.x,Cug.y,Cug.z,titl,ibar,iaoff)
          else
           plot(xvero,Etcf)
           title(titl)
          end
         end
        else
        subplot(nsub,nsub2,nplo+(nsub2-1)*(nplo-1))
         if iFF==1
          if nsub2==2
           if iFFte==0
            titl=[' FF at ',num2str(z0c),'cm; scales in cm '];
           else
            titl=[' FF radiation pattern '];
           end
           if i2D==3 & suee>0
            surf(X,Y,(Ef).^2),

% x0=-max(max(X)); x1=-x0; y0=5; y1=5.5;
% fi=find(( X>x0&X<x1) & (X>x0&X<x1) );
% Ef2=Ef.^2;
% [du,im]=max(max(Ef2));
% fim=find(Ef2==du)
% fim0=fim(1);
% [X(fim0) Y(fim0)]
%  roM=max(roxt)*1e-4;
%  Xi=linspace(-1,1,200)*roM;
%  Yi=z0c*tan(X(fim0)*pi/180);
% Fi=atan(Xi/Yi);
% fim=find(Fi<0);
% Fi(fim)=Fi(fim)+2*pi;
% Ri=sqrt(Yi.^2+Xi.^2);
% rX=roxt*1e-4;
% figure, polar(Fi,Ri)
%
% Z=interpn(rX,fian,log(Ef2),Ri,Fi,'spline');
% FF=exp(Z);
% Xid=atan(Xi/z0c)*180/pi;
% figure, plot(Xid,FF,'r'), hold on, plot(X(fi),Ef2(fi),'.')
%

% figure, plot(Y(fi),Ef(fi).^2,'.')

%           figure
%            surf(X,Y,abs(Hfy).^2),
%           figure
%            surf(X,Y,abs(Efx).^2),
%           figure
%            surf(X,Y,abs(Ef).^2),
            shading('interp'), view(0,90),
            axis square, axis equal, grid,
            axis([-1 1 -1 1]*axlim/2),
            title(titl)
%         if ieftot==0
%          Eft=Ef.^2;
%         elseif ieftot==1
%          Eft=abs(Efx).^2;
%         else
%          Eft=abs(Efy).^2;
%         end         
 tet=X(:,1)/180*pi;
 dtet=tet.*[0; diff(tet)];
         EF{imazver,nso}=Ef.^2;
         EFFy{imazver,nso}=abs(Efy).^2;
%         EF(:,:,imazver)=Efx;
            if  imazver==1
             effma=(max(max(Ef.^2)));
%             dted=diff(X(:,1));
%             fate=pi/180;
%             dte=[dted(1); dted].*X(:,1)*fate^2; 
%             dfi=ones(size(fian))*diff(fian(1:2));
             Ecaf=mean(abs(Ef).^2,2);
             Enf=dtet'*Ecaf;
%             keyboard
%             effmaint=dte'*Ef.^2*dfi';
             effmaint=Enf;
    
             effintv(imazver,nso)=Enf;
            else
%             effintv(imazver)=dte'*Ef.^2*dfi';
            
             Ecaf=mean(abs(Ef).^2,2);
             Enf=dtet.'*Ecaf;
             effintv(imazver,nso)=Enf;
%             effmaint=sum(sum(Ef.^2));
            end
            if imazver==1000
             figure(ha)
             subplot(3,2,imazver)
%             surf(X,Y,(Ef).^2/effma), colorbar
             surf(X,Y,Efx/sqrt(effma)), colorbar
                         shading('interp'), view(0,90),
	                 axis square, axis equal, grid,
	                 axis([-1 1 -1 1]*axlim/2),
	                 if imazver==1
	                  tit=[' Totale; Integrale =1 '];
	                  tit=[' Totale'];
	                 else
                         effint=dte'*Ef.^2*dfi';	                 
	                  tit=[' Contributi ', num2str(Mazv{imazver}),'; Integrale =',num2str(effint/effmaint,3)];
	                  tit=[' Contributi ', num2str(Mazv{imazver}),];
	                 end
	                 title(tit)
	                 pausak
            end
                        
            
           end


          end
         end

        end

       else


         if iFF==1
          if nsub2==2
           if iFFte==0
            titl=[' FF at ',num2str(z0c),'cm; scales in cm '];
           else
            titl=[' FF radiation pattern '];
           end
           if i2D==3
            surf(X,Y,(Ef).^2),
            shading('interp'), view(0,90),
            axis square, axis equal, grid,
            axis([-1 1 -1 1]*axlim/2),
            title(titl)
           end
          end
         end


       end

%       keyboard

    ifp=ifps;
    end

   else
    graph
    if ifp==-10 & iLP==0
%     pol_rpo
    end
   end
   % fine plot campi
    tet=X(:,1)/180*pi;
    dtet=tet.*[0; diff(tet)];
            EF{imazver,nso}=Ef.^2;
            EFFy{imazver,nso}=abs(Efy).^2;
   %         EF(:,:,imazver)=Efx;
               if  imazver==1
                effma=(max(max(Ef.^2)));
   %             dted=diff(X(:,1));
   %             fate=pi/180;
   %             dte=[dted(1); dted].*X(:,1)*fate^2; 
   %             dfi=ones(size(fian))*diff(fian(1:2));
                Ecaf=mean(abs(Ef).^2,2);
                Enf=dtet'*Ecaf;
   %             keyboard
   %             effmaint=dte'*Ef.^2*dfi';
                effmaint=Enf;
       
                effintv(imazver,nso)=Enf;
               else
   %             effintv(imazver)=dte'*Ef.^2*dfi';
               
                Ecaf=mean(abs(Ef).^2,2);
                Enf=dtet.'*Ecaf;
                effintv(imazver,nso)=Enf;
   %             effmaint=sum(sum(Ef.^2));
            end
   
%   'VERT ', keyboard
   if SPP==10 & ifp==-10
     fase_arm
     figure(fgsav)
   end
   
%       if ifp==-10 & M2<M2_max, pausak, end
       if ifp==-10, pausak, end


 rg=rr;
 if idyn>=1
  Ef_sav=Ef;
  if i2D==3
   Ec_at=Etb(:,:,2);
  else
   Ec_at=Etb(:,2);
  end
  Ec_as=Ec_at;
  Int_cb=xdx*Ec_as*defi;

  if N==0
%  if iplan==1
   Ec_atc=Ec_at;
   Ec_atc(fiGa)=0;
  else
   Nga=[N; zeros(length(x)-length(ro_in),1)];
   Nmat=repmat(Nga,1,length(fian));
%   NNor=xdx*Nga*2*pi;
%   Are=pi*avero^2;
   NNor=1;
   Are=1;
   Ec_atc=Ec_at.*Nmat/NNor*Are;
  end

  Int_at=xdx*Ec_atc*defi;
  Gam_at=Int_at/Int_cb;

%' gam', keyboard



  cmi=c*1e6;
  omega=cmi*k0*(1+fris);
  Perdvol=-omega*sum(Perd)/confz;
  if iplan==0
   Glostotal=2*gg0*Gam_at+Perdvol*(1-iPERD_vol);
  else
   Glostotal=2*gg0+Perdvol*(1-iPERD_vol);
  end
  

  
  Lsov(imod)=Glostotal*confz;
%  keyboard
%  keyboard
%  keyboard
%  keyboard
%  keyboard

  Ec_su=Etcsum;
  Int_su=xdx*Ec_su*defi;

  fa4=Int_su/Int_cb*fatqw;
  Gam_v(imod)=Gam_at;

  if i2D==3
   Ec_at=Etf(:,:,2);
  else
   Ec_at=Etf(:,2);
  end
%  Int_cb=xdx*Ec_as*defi;
  
  Int_cf=xdx*Ec_at*defi;

%   disp(' Completamente planare '), keyboard
  if ifp==-10
%   keyboard
   pausak
  end

  if ifp>-4 | ifp==-10
   disp('  ')
   disp('  ')
   disp(' Completamente planare ')
   disp('  ')
  end
  rg=rr;
  rg=nmean;
  CP=c/rg/(2*Lef*confz);
  CP0=c/rg/(2*Lef);
  Gsup=fTras*CP
  Ginf=fTrasinf*CP
  Gper=(fPinf+fPsup)*CP/2
  Gplanar=Gsup+Ginf+Gper

  if ifp>-4 | ifp==-10
   disp('  ')
   disp('  ')
   disp('  ')
   disp(' Planare corretto  da integrali')
   disp('  ')
  end

  if i2D==3
   Ef=Etf(:,:,1);
   Eb=Etb(:,:,1);
   Po=Poi(:,:,1);
  else
   Ef=Etf(:,1);
   Eb=Etb(:,1);
   Po=Poi(:,1);
  end
  Euf=Ef;
  Eub=Eb;
  Pou=Po;

  Int_f=xdx*Ef*defi;
  Int_b=xdx*Eb*defi;
  Int_P=xdx*Po*defi;
  Int_fsu=Int_f;
  Int_bsu=Int_b;
 for ipo=1:icsmax
  Int_Po(ipo)=xdx*Poi(:,:,ipo)*defi;
 end
 if icsmax==5
  IP_ord=real(Int_Po([1 4 2 5 3]));
 else
  %IP_ord=real(Int_Po([1 2 5 3]));
  IP_ord=real(Int_Po([1 2 4 3]));
 end 
 dP=[diff([0 IP_ord 0])];
 [maP,imadu]=max(abs(dP));
 dP(imadu)=0;
 dP=dP(find(dP~=0))/maP;
 dP_perc=dP*100;
 dPv=dP*glostotal;
 fiesc=find(abs(dPv)<.002*max(dPv));
 dPv(fiesc)=0;
 %'Poi', keyboard
 
  Goutsup=(Int_b-Int_f)*CP/Int_cb
  GoutsupP=Int_P*CP/Int_cb
%  Goutsup=(Int_b+Int_f)*CP/Int_cb
%   'Perdsup ', keyboard

    kmo=2;
    Andu=Afz(:,kmo)+Afzf(:,kmo);
    cam_val;
    Eanuo=Etot;
    Int_at=xdx*Eanuo*defi;
    Teff=(Int_b-Int_f)/Int_at;



  if i2D==3
   Ef=Etf(:,:,3);
   Eb=Etb(:,:,3);
  else
   Ef=Etf(:,3);
   Eb=Etb(:,3);
  end
  Int_f=xdx*Ef*defi;
  Int_b=xdx*Eb*defi;
  Int_fso=Int_f;
  Int_bso=Int_b;  
  Ebb=Eb; 
  Ebf=Ef; 
  Goutinf=(Int_f-Int_b)*CP/Int_cf
%  Goutinf=(Int_f+Int_b)*CP/Int_cf

  if istmet==0
   istmet1=[];
  else
   istmet1=istmet;
  end
  if length(istmet1)==1
   if i2D==3 
    Ef=Etf(:,:,istmet);
    Eb=Etb(:,:,istmet);
    Po=Poi(:,:,istmet);
   else
    Ef=Etf(:,istmet);
    Eb=Etb(:,istmet);
    Po=Poi(:,istmet);
   end
   
   Int_f=xdx*Ef*defi;
   Int_b=xdx*Eb*defi;
   Int_P=xdx*Po*defi;
   Pom=Po;
  Emb=Eb; 
  Emf=Ef; 
  Int_fme=Int_f;
  Int_bme=Int_b;   
  
 % Goutsup=(Int_b-Int_f)*CP/Int_cb 
  
  Gtransmet=(Int_b-Int_f)*CP/Int_cb;
  GtransmetP=Int_P*CP/Int_cb;
  
   Perdring=Gtransmet-Goutsup
%   'Perdring ', keyboard
  else
   Perdring=0;
  end

  Perdrad=Goutsup+Goutinf;
    Perdrad_s=Goutsup;
    Perdrad_i=Goutinf;

  Perdvol_u=-omega*Perd(1)/confz;
  Perdvol_d=-omega*Perd(2)/confz;
  Gapprox=(Perdvol*iPERD_vol+Perdrad+Perdring);
   disp(' confronto approx  numerico ')
   [Gapprox glosout*Gam_at]

  if ifp==-10
%   keyboard
   pausak
  end

  if ifp>-4 | ifp==-10
   disp('  ')
   disp('  ')
   disp('  ')
   disp('[total= pred_vol + out_sup + out_inf]');
   disp('  ')
   disp(' Planar ');
   [Gplanar Gper Gsup Ginf ]

   disp('  ')
   disp(' integral method ');
   [Gapprox  Perdvol Goutsup Goutinf]

   disp(' Confronto completamente planare: verifica '),
  end


  Wv(imod)=W;
  czv(imod)=confztot;
  fqwv(imod)=fatqw;
  fa4v(imod)=fa4;
  Lefv(imod)=Lef;

  Tefv(imod)=fTras;
  Tefiv(imod)=fTrasinf;
  gu=fTras*CP0;
  gi=fTrasinf*CP0;



%  taut(imod)=1/(glosout*Gam_at*confz);
%  tauu(imod)=1/((Goutsup+Perdvol_u)*CP0/CP);
%  taub(imod)=1/((Goutinf+Perdvol_d)*CP0/CP);
  taut(imod)=1/Lsov(imod);
%  tauu(imod)=1/((Goutsup)*CP0/CP);
%  taub(imod)=1/((Goutinf)*CP0/CP);

  teff(imod)=Teff;

%  al_c=al_cav(fou(nso,:),aou(nso,:),gou(nso,:),ze,gso);
  alca(imod)=0;

%  tauu(imod)=1/((Goutsup)*confz);
%  taub(imod)=1/((Goutinf)*confz);
%  pvol(imod)=1/(Perdvol*confz);
%  pmet(imod)=1/(Perdring*confz);
  confz=uL(1);
    tauu(imod)=1/(dPv(1)*confz);
    taub(imod)=1/(dPv(end)*confz);
    pvol(imod)=1/(sum(dPv(end-2:end-1))*confz);
    if icsmax==5
     pmet(imod)=1/(dPv(2)*confz);
    else
     pmet(imod)=0;
    end
  
  dP_old=[Goutsup Perdring Perdvol_u Perdvol_d Goutinf];
  
 if ifp==-10
  pux=0;
  if icsmax==5
   ipux=istfi.t;
   pux=[pux sum(Litot(ipux))];
  end
  fiqw=find(iauto(:,1)==2)-2;
  Lprod=Litot.*fmlstot(:,2);
  Lqw=sum(Lprod(1:fiqw));
  Lto=sum(Lprod);
  pux=[pux Lqw+[0 Litot(fiqw+1)] Lto]; 
  
  figure, plot(pux,IP_ord'/maP,'o'),
  xlabel(' long coord (micron)')
  ylabel(' normalized loss distribution ')
  ' ver perdite', keyboard
 end
  global irigau

  if imod==1
   if length(irigau)>0
    if irigau==1
     save riga
%    ' iriga ', keyboard
%    rigasub
    end
   end
  end

%  disp(' Controllo S21  mio ')
%  keyboard

  if ifp>-4 | ifp==-10
   disp('  ')
   disp('  ')
   disp('  ')
   disp('[total: numerico, integrale, planare]');
   disp('  ')
   [glostotal*Gam_at Gapprox Gplanar]

   disp(' camfull: verifica prova perdite ')
  end

   itappo=0;   %calcolo varie perdite con eig
   if itappo==1
      tappo
   end  %itappo
  if  ifp==-10
%   keyboard
   pausak
  end
  Ef=Ef_sav;
 end  %idyn


else   %iskim
    imod=imod-1;
end

if iskim==0
 if ifp==-10
%  pausak
 end
% if ifp==-4
%  raf=[nomeFs(1:end-4),num2str(iLP),num2str(mm)];
%  indi=num2str(fig_ind.a);
%  if length(fig_ind.b)>0
%   indi=[indi num2str(fig_ind.b)];
%  end
%  indi=[indi num2str(imod)];
%  stri=[' fig.ind = ',indi];
%  subplot(nsub,nsub2,1)
%  text(-.3,1.2+shi,stri);
%  eval(['hgsave(',num2str(fgsav),',''',raf,indi,''');']);
%  close all
% end

%  Plot.Ap{imod,imm}=Anout;
%  Plot.ApQ{imod,imm}=AnQW;
%  Plot.parmod{imod,imm}=[glosout lambda/(1+ze) ze pola alim nK_dis];
%  Plot.XP{imod,imm}=XP;
%  Plot.YP{imod,imm}=YP;
%  Plot.X{imod,imm}=X;
%  Plot.Y{imod,imm}=Y;
%  Plot.E2xo{imod,imm}=E2xo;
%  Plot.E2xp{imod,imm}=E2xp;
%  if iLP==0
%   Plot.E2yo{imod,imm}=E2yo;
%   Plot.E2yp{imod,imm}=E2yp;
%   Plot.E2zo{imod,imm}=E2zo;
%  end
%  Plot.Ef{imod,imm}=Ef;
%  Plot.Cug{imod,imm}=Cug;
 if imazver==1
  Plot.Ef{imod}=Ef;
  Plot.Efx{imod}=Efx;
  Plot.Efy{imod}=Efy;
  Plot.Cug{imod}=Cug;

  Plot.Ap{imod}=Anout;
  Plot.ApQ{imod}=AnQW;
  Plot.parmod{imod}=[glosout lambda/(1+ze) ze pola alim nK_dis];
  Plot.XP{imod}=XP;
  Plot.YP{imod}=YP;
  Plot.X{imod}=X;
  Plot.Y{imod}=Y;
  Plot.E2xo{imod}=E2xo;
  Plot.E2xp{imod}=E2xp;
  global ikr
   FF.ikr=ikr;
   FF.segem=segem;
      FF.mbv=mbv;
   FF.numodi=numodi;
   FF.numodiacc=numodiacc;
   FF.Pus=Pus;
   FF.iredmat=iredmat;
   FF.Anout=Anout;
   FF.Gas=Gas;
   FF.Gad=Gad;
   FF.iFFsect=iFFsect;
   FF.lambda=lambda;
   FF.ze=ze;
   FF.nv=nv;
   FF.rr=rr;
   FF.KKv=KKv;
   FF.KKt=KKt;
   FF.besm=besm;
   FF.Mvefm0=Mvefm0;
   if iLP==0
    FF.besp=besp;
    FF.besz=besz;
    FF.Mvez0=Mvez0;
    FF.Mvhz0=Mvhz0;
    FF.Mvefp0=Mvefp0;
    FF.Mvegm0=Mvegm0;
    FF.Mvegp0=Mvegp0;
   end
   Znor=Zve;
   FF.Znor=Znor;
   FF.pes=pes;
   FF.iLP=iLP;
   FF.Pus=Pus;
   FF.iLP=iLP;
   FF.xvero=xvero;
   FF.aiat=aiat;
   FF.Litot=Litot;
   FF.aitot=aitot;
   FF.nitot=nitot;
   FF.iauto=iauto;
   FF.lambda=lambda;
   FF.fmlstot=fmlstot;
   FF.Pusas=Pusas;
   FF.XP=XP;
   FF.YP=YP;
   FF.ifp=ifp;

   global iphase_front isav_nfcut  Ps


%'Az', keyboard
   fi_save=fi;
   
   if abs(isav_Az)>=1 & imod<=3
%   ifp=-11
    if ifp==-10
    ' sono qui prima di sub_Ez'
    keyboard
    end
%   ' sono qui prima di sub_Ez', keyboard

%ipoS=ipolar;
%ipolar=pola;
%    iLP=iLPr;
    polasa=pola;
    polsav=pol;
    sub_Ez
    pola=polasa;
    pol=polsav;
    FF.Nz=Nz;
%   else
%    ' Acoz'
%   keyboard
%    FF.Acoz=
   end


   if isav_nfcut==1
    global lp1 lp2 lp3
    Nofz=[nomeFe,'_',num2str(lp1),num2str(lp2),'_',num2str(imod)];
    eval(['save ',Nofz])
    [' salvato file:  --->    ',Nofz], keyboard
   end

   if length(iphase_front)==0
    iphase_front=0;
   end

   if iphase_front==1
    H0=1.5;
    [D0,R0,fav,faf]=ffl_sub(FF,H0);
    FF.D0=D0;
    FF.R0=R0;
    FF.fav=fav;
    FF.faf=faf;
    if ifp~=-4
     figure, plot(fav(:,1),fav(:,2),'r.',faf(:,1),faf(:,2))
    ' ff_sub'
     keyboard
    end
   end

   Plot.FF{imod}=FF;
   if iLP==0
      I=abs(abs(E2xp.^2)+abs(E2yp.^2));
   else
      I=abs(abs(E2xp.^2));
   end
%   Plot.E2xp{imod}=I;

  if iLP==0
   Plot.E2yo{imod}=E2yo;
   Plot.E2yp{imod}=E2yp;
   Plot.E2zo{imod}=E2zp;

   I=abs(abs(Efx.^2)+abs(Efy.^2));
   V=j*(Efx.*conj(Efy)-conj(Efx).*Efy);
   In=max(max(I));
   Dci=V./In;
%   Plot.Dc{imod}=Dci;
   Plot.E2zp{imod}=Dci;
  end




%  Plot.Ef{imod}=Ef;
%  Plot.Cug{imod}=Cug;
  ifouu=0;
if ifouu==0  
%  if igouadd==1
  Plot.gou{imod}=gou;
  Plot.aou{imod}=aou;
  Plot.fou{imod}=fou+dla_fre;
%  igouadd=0;  
%  end
else  
  Plot.gou{imod}=ga;
  Plot.aou{imod}=ya1;
   fi=fi_save;
   fa1d=fou(nso,ps);
   fa1d=fa1d(fi);
  Plot.fou{imod}=fa1d;
end  
  
  Plot.ze{imod}=ze+dla_fre;
  Plot.gg0{imod}=gg0;

%  ' Plot ', keyboard

%  KK numodi pasnu lbv kcav0 x Nx mbv KKt fian
end  %imazver

end

global iffcut fun_ff


if length(iffcut)>0
 if iffcut==1 & imod==1
  if length(fun_ff)==0
   ff_cut
  else
   eval(fun_ff)
   if ifp~=-4
    pausak
   end
  end
%  cut_e
 end
end

iLP=iLPr;