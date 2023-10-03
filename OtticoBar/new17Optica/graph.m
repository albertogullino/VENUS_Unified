
%aax=0.8*max(max(XP));
if ifp<=-10
    ifps=ifp;
    ifp=1;
    fgsav=figure;

%     if ~exist('ipolari')
      if iLP==0
       nsub=3
       pos=20;
      else
       nsub=2
       pos=20;
      end
%     else
%       nsub=3
%     end
%     pos
%'qui'
%keyboard
 if pola==1
  pograp(1)=pograp(1)+pos*sinc;
  set(gcf,'Position',pograp);
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
  pogram(1)=pogram(1)+pos*sinc;
  set(gcf,'Position',pogram);
 end


%      if pola==1
%       pograp(1)=pograp(1)+pos*sinc;
%       set(gcf,'Position',pograp);
%      else
%       pogram(1)=pogram(1)+pos*sinc;
%       set(gcf,'Position',pogram);
%      end


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
      text(fix(length(Ap)/2.5),1.27+shi,stri);
      stri=[' wav = ',num2str(lambda/(1+ze),'%0.4e')];
      text(fix(length(Ap)/2.5),1.17+shi,stri);
      stri=[' freq = ',num2str(ze,'%0.4e')];
      text(fix(length(Ap)/2.5),1.07+shi,stri);


       titl=' E_x  Output';
       nplo=2;
       subplot(nsub,nsub2,nplo+(nsub2-1)*(nplo-1)+1)
%       subplot(nsub,nsub2,3)
       ibar=1;

       if i2D==3
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
          map_fnew(XP,YP,E2mp,aax,Cug.x,Cug.y,Cug.z,titl,ibar,iaoff)
          iaoff=1;
         end
        end

       end

        subplot(nsub,nsub2,2)
       if iEz==1
        titl=' E_z  Output';
        map_fnew(XP,YP,E2zo,aax,Cug.x,Cug.y,Cug.z,titl,ibar,iaoff)
       else
         if iFF==1
          if nsub2==2
           if iFFte==0
            titl=[' FF at ',num2str(z0c),'cm; scales in cm '];
           else
            if ~exist('parlop')
             titl=[' FF radiation pattern '];
            else
             titl=[' FF radiation pattern ',num2str(parlop)];
            end
           end
           if i2D==3
            surf(X,Y,abs(Ef).^2),
            shading('interp'), view(0,90),
            axis square, axis equal, grid,
            axis([-1 1 -1 1]*axli/2),
            title(titl)
           end
          end
         end


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
            if ~exist('parlop')
             titl=[' FF radiation pattern '];
            else
             titl=[' FF radiation pattern ',num2str(parlop)];
            end
           end
           if i2D==3
            surf(X,Y,abs(Ef).^2),
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
            if ~exist('parlop')
             titl=[' FF radiation pattern '];
            else
             titl=[' FF radiation pattern ',num2str(parlop)];
            end
           end
           if i2D==3
            surf(X,Y,abs(Ef).^2),
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




%        figure
%        map_fnew(XP,YP,E2xo,aax,Cug.x,Cug.y,Cug.z,titl,ibar,iaoff)
%        figure
%        map_fnew(XP,YP,E2xp,aax,Cug.x,Cug.y,Cug.z,titl,ibar,iaoff)
