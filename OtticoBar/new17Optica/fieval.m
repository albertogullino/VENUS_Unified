
iEz=1;
fatEz=rr^2;
if iLP==1
 iEz=0;
end
'Ez'
keyboard
imod=imod+1
iaoff=1;
iskim=0;        % serve per LP e modo fondamentale tipo sin(0*fi)
icut=0;
if ~exist('icaplo')
 icaplo=1;
end
if length(icaplo)==0
 icaplo=1;
end

  glosout=2*gg0;
  glostotal=2*gg0;
  Gsimmi=gg0;
  fris=ze;
  freq=fris;
  ADoutis=An0;
  ADouti=An0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% decido quali punti prendere dell'autovettore (puAc)

    taglio

% FINE decido quali punti prendere dell'autovettore
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calcolo tutti i campi nelle sezioni che servono


 calc_new


%  FINE calcolo tutti i campi nelle sezioni che servono
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% assegno i vari parametri modali (se il modo e` buono: iskim=0)


if iskim==0

   glostotal=glosout;

   ADc(:,imod)=Afz(:,2);
   ADout(:,imod)=Afz(:,1);
   Gsov(imod)=glostotal;
   Fsov(imod)=ze;
   Nazim(imod)=nuazi;
   Nrad(imod)=mrad;
   tyE(imod)=1;

   if idyn>=1
    TL_ef
    W=0;
   end

   if i2D==3

     s=size(Etot);
     p1=1:s(1);
     p2=1:s(2);
     EDc.x(p1,p2,imod)=Etx(:,:,2);
     EDout.x(p1,p2,imod)=Etx(:,:,1);
     E2xo=Etx(:,:,1);
     E2xp=Etx(:,:,icaplo);

     if iLP==0
      EDc.y(p1,p2,imod)=Ety(:,:,2);
      EDout.y(p1,p2,imod)=Ety(:,:,1);
      E2yo=Ety(:,:,1);
      E2zo=Etz(:,:,1)*fatEz;
      E2yp=Ety(:,:,icaplo);
     end

     if iLP==0
      if max(max(E2y))>max(max(E2x))
       tyE(imod)=-1;
      end
     end

   else

     EDc.x(:,imod)=Etx(:,2);
     EDout.x(:,imod)=Etx(:,1);
     if iLP==0
      EDc.y(:,imod)=Ety(:,2);
      EDout.y(:,imod)=Ety(:,1);
     end

   end

      Ap=abs(Anout)/max(abs(Anout));
      ApQ=abs(AnQW)/max(abs(AnQW));

%    if ifp<=-10 | ifp==-4
    if ifp<=-10
    ifps=ifp;
    ifp=1;
     fgsav=figure;

      if iLP==0
       nsub=3
       pos=40;
      else
       nsub=2
       pos=10;
      end
      if pola==1
       set(gcf,'Position',pograp);
       pograp(1)=pograp(1)+pos*sinc;
       if sinc<1
        if pograp(1)<50
         pograp(1)=950;
         pograp(2)=100;
        end
       else
       if pograp(1)>950
        pograp(1)=50;
        pograp(2)=100;
       end
       end
      else
       set(gcf,'Position',pogram);
       pogram(1)=pogram(1)+pos*sinc;
      end



      subplot(nsub,nsub2,1)



      if pola==1
       plot(Ap,'r.-'), hold on, plot(ApQ,'w') ;
      else
       plot(Ap,'g.-'), hold on, plot(ApQ,'w') ;
      end
      if iLP==0
       shi=0;
      else
       shi=-.2;
      end
      stri=[' gain = ',num2str(glosout,'%0.4e')];
      text(fix(length(Anout)/2.5),1.27+shi,stri);
      stri=[' wav = ',num2str(lambda/(1+ze),'%0.4e')];
      text(fix(length(Anout)/2.5),1.17+shi,stri);
      stri=[' freq = ',num2str(ze,'%0.4e')];
      text(fix(length(Anout)/2.5),1.07+shi,stri);


       titl=' E_x  Output';
       nplo=2;
       subplot(nsub,nsub2,nplo+(nsub2-1)*(nplo-1)+1)
%       subplot(nsub,nsub2,3)
       ibar=1;

       if i2D==3
        aax=0.8*max(xvero);
        map_fnew(XP,YP,E2xo,aax,Cug.x,Cug.y,Cug.z,titl,ibar,iaoff)
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
         titl=' E_x  QW';
         nplo=2;
         subplot(nsub,nsub2,nplo+(nsub2-1)*(nplo-1))
         ibar=1;
         iaoff=0;
         if i2D==3
          aax=0.8*max(xvero);
          map_fnew(XP,YP,E2xp,aax,Cug.x,Cug.y,Cug.z,titl,ibar,iaoff)
          iaoff=1;
         else
          plot(xvero,Etcf)
          title(titl)
         end
        end
       else
        if nsub2==2
         titl=titlm;
         nplo=2;
         subplot(nsub,nsub2,nplo+(nsub2-1)*(nplo-1))
         ibar=1;
         iaoff=0;
         if i2D==3
          aax=0.8*max(xvero);
          map_fnew(XP,YP,E2mp,aax,Cug.x,Cug.y,Cug.z,titl,ibar,iaoff)
          iaoff=1;
         end
        end

       end

       subplot(nsub,nsub2,2)
       if iEz==1
        titl=' E_z  Output';
        aax=0.8*max(xvero);
        map_fnew(XP,YP,E2zo,aax,Cug.x,Cug.y,Cug.z,titl,ibar,iaoff)
       end


       if iLP==0 & i2D==3
%        titl=' |E|^2 ';
%        subplot(nsub,1,3)
        titl=' E_y Output';
        nplo=3;
        subplot(nsub,nsub2,nplo+(nsub2-1)*(nplo-1)+1)
%        subplot(nsub,nsub2,5)
         map_fnew(XP,YP,E2yo,aax,Cug.x,Cug.y,Cug.z,titl,ibar,iaoff)


        if iEz==0
        subplot(nsub,nsub2,nplo+(nsub2-1)*(nplo-1))
         if nsub2==2
          titl=' E_y  QW';
          nplo=3;
          ibar=1;
          if i2D==3
           aax=0.8*max(xvero);
           map_fnew(XP,YP,E2yp,aax,Cug.x,Cug.y,Cug.z,titl,ibar,iaoff)
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
           if i2D==3
            surf(X,Y,(Ef).^2),
            shading('interp'), view(0,90),
            axis square, axis equal, grid,
            axis([-1 1 -1 1]*axli/2),
            title(titl)
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
            axis([-1 1 -1 1]*axli/2),
            title(titl)
           end
          end
         end


       end

       if ifp==-10, pausak, end
%       keyboard

    ifp=ifps;
    end

 if idyn>=1

  if i2D==3
   Ec_at=Etb(:,:,2);
  else
   Ec_at=Etb(:,2);
  end

  Int_c=xdx*Ec_at*defi;
  Ec_at(fiGa)=0;
  Int_at=xdx*Ec_at*defi;
  Gam_at=Int_at/Int_c;

  Ec_su=Etcsum;
  Int_su=xdx*Ec_su*defi;

  fa4=Int_su/Int_at*fatqw;
  Gam_v(imod)=Gam_at;



  if ifp>-4
   disp('  ')
   disp('  ')
   disp(' Completamente planare ')
   disp('  ')
  end

  CP=c/rr/(2*Lef*confz);
  Gsup=fTras*CP
  Ginf=fTrasinf*CP
  Gper=(fPinf+fPsup)*CP/2
  Gplanar=Gsup+Ginf+Gper

  if ifp>-4
   disp('  ')
   disp('  ')
   disp('  ')
   disp(' Planare corretto  da integrali')
   disp('  ')
  end

  if i2D==3
   Ef=Etf(:,:,1);
   Eb=Etb(:,:,1);
  else
   Ef=Etf(:,1);
   Eb=Etb(:,1);
  end

  Int_f=xdx*Ef*defi;
  Int_b=xdx*Eb*defi;

  Goutsup=-(Int_f-Int_b)*CP/Int_c
  
  Perdrad_s=Goutsup;


  if i2D==3
   Ef=Etf(:,:,3);
   Eb=Etb(:,:,3);
  else
   Ef=Etf(:,3);
   Eb=Etb(:,3);
  end
  Int_f=xdx*Ef*defi;
  Int_b=xdx*Eb*defi;
  Goutinf=(Int_f-Int_b)*CP/Int_c
  
  Perdrad_i=Goutsup;


  if istmet==0
   istmet1=[];
  else
   istmet1=istmet;
  end
  if length(istmet1)==1
   if i2D==3
    Ef=Etf(:,:,istmet);
    Eb=Etb(:,:,istmet);
   else
    Ef=Etf(:,istmet);
    Eb=Etb(:,istmet);
   end
   Int_f=xdx*Ef*defi;
   Int_b=xdx*Eb*defi;
   Gtransmet=-(Int_f-Int_b)*CP/Int_c;
   Perdring=Gtransmet-Goutsup
  else
   Perdring=0;
  end

  Perdrad=Goutsup+Goutinf;
  cmi=c*1e6;
  omega=cmi*k0*(1+freq);
  Perdvol=-omega*Perd/confz
  Gapprox=(Perdvol+Perdrad+Perdring);
   disp(' confronto tot / approx ')
   [Gapprox glosout]
   keyboard

  if ifp>-4
   disp('  ')
   disp('  ')
   disp('  ')
   disp('[total pred_vol out_sup  out_inf]');
   disp('  ')
   disp(' Planar ');
   [Gplanar Gper Gsup Ginf ]

   disp('  ')
   disp(' integral method ');
   [Gapprox  Perdvol Goutsup Goutinf]

   disp(' Confronto completamente planare: verifica '),
  end


  Wv(imod)=W;
  czv(imod)=confz;
  fqwv(imod)=fatqw;
  fa4v(imod)=fa4;
  Lefv(imod)=Lef;
  Tefv(imod)=fTras;
  Tefiv(imod)=fTrasinf;
  gu=fTras*c/rr/(2*Lef);
  gi=fTrasinf*c/rr/(2*Lef);
%  gtot(imod)=glostotal*confz;

  gtot(imod)=Goutsup/CP;

  Perdu(imod)=fPsup*CP/2;
  Perdb(imod)=fPinf*CP/2;
  Perdb(imod)=fPinf*CP/2;
  disp(' Controllo S21  mio ')
  keyboard

  if ifp>-4
   disp('  ')
   disp('  ')
   disp('  ')
   disp('[total: numerico, integrale, planare]');
   disp('  ')
   [glostotal Gapprox Gplanar]

   disp(' camfull: verifica prova perdite ')
  end

   itappo=0;   %calcolo varie perdite con eig
   if itappo==1
      tappo
   end  %itappo
keyboard
 end  %idyn


else   %iskim
    imod=imod-1;
end

if iskim==0
 if ifp==-10
  pausak
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

  Plot.Ap{imod,imm}=Anout;
  Plot.ApQ{imod,imm}=AnQW;
  Plot.parmod{imod,imm}=[glosout lambda/(1+ze) ze];
  Plot.XP{imod,imm}=XP;
  Plot.YP{imod,imm}=YP;
  Plot.X{imod,imm}=X;
  Plot.Y{imod,imm}=Y;
  Plot.E2xo{imod,imm}=E2xo;
  Plot.E2xp{imod,imm}=E2xp;
  if iLP==0
   Plot.E2yo{imod,imm}=E2yo;
   Plot.E2yp{imod,imm}=E2yp;
   Plot.E2zo{imod,imm}=E2zo;
  end
  Plot.Ef{imod,imm}=Ef;
  Plot.Cug{imod,imm}=Cug;



%  KK numodi pasnu lbv kcav0 x Nx mbv KKt fian

end
