function map_fiel(XP,YP,E,aax,xcui,ycui,zcui,titl,ibar,iaoff)
col(1,1)='w';
col(1,2)='k';
col(1,3)='y';
%'in map'
%keyboard
Emoto=real(E);

    surf(XP,YP,Emoto),
    shading('interp'),
    view(0,90),
    colormap('gray')
    axis square,
    axis equal,
    axis([-1 1 -1 1]*aax),
    if iaoff==1
     axis off
    end
    fczu=1;
    title(titl),
    if abs(ibar)==1
     if ibar==-1
     colorbar('horiz')
     else
     colorbar
     end
    end
    si=size(xcui);
   if prod(si)>1
    hold on,
    for jj=1:si(2)
     shav=xcui(1,jj);
     npo=xcui(2,jj);
     pux=3:2+npo;
     xcu=xcui(pux,jj);
     ycu=ycui(pux,jj);
     zcu=zcui(pux,jj);
     if shav~=5
      ph=plot3(xcu,ycu,zcu*fczu,'w');
      set(ph,'linewidth',1)
      if shav==6
       set(ph,'linewidth',1)
      end
     else
      lab=xcui(3+npo,jj);
      if lab==1
       ph=plot3(xcu,ycu,zcu*fczu,'w');
      elseif lab==2
       ph=plot3(xcu,ycu,zcu*fczu,'k');
      elseif lab==3
       ph=plot3(xcu,ycu,zcu*fczu,'y');
      end
      set(ph,'linewidth',1)
     end
    end
   end
