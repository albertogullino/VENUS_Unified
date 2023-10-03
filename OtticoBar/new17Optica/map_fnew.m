%function map_fiel(XP,YP,E,aax,xcui,ycui,zcui,titl,ibar,iaoff,icol,ilini,itexti)
%

function map_fiel(XP,YP,E,aax,xcui,ycui,zcui,titl,ibar,iaoff,icol,ilini,itexti)

%'entro vor', keyboard
col(1,1)='w';
col(1,2)='k';
col(1,3)='y';
Emoto=real(E);
warning off
if nargin<=11
 iaoff=1;
end

    surf(XP,YP,Emoto),
%    keyboard
%    map(Emoto,XP,YP),
    shading('interp'),
if nargin>=11
 if icol==0
  colormap('gray')
 end
end

ilin=1.5;
if nargin>=12
 if length(ilini)>0
  ilin=ilini;
 end
end
itext=0;
if nargin>=13
 itext=itexti;
end
    view(0,90),
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
     h=colorbar('horiz');
     else
     h=colorbar;
     end
     if itext>0
      set(h,'fontsize',itext)
     end
    end
%    'mat fnew', keyboard
    if ilin~=0
    si=size(xcui);
   if prod(si)>1
%    if si(2)>4
%     si2=[1 si(2)];
%    else
     si2=[1:si(2)];
%    end
    hold on,
    for jj=si2
     shav=xcui(1,jj);
     npo=xcui(2,jj);
     pux=3:2+npo;
     xcu=xcui(pux,jj);
     ycu=ycui(pux,jj);
     zcu=zcui(pux,jj);
     %if shav==1
      co='w';
     %else
     % co='g';
     %end 
     if shav~=5
%      if shav~=6
      ph=plot3(xcu,ycu,zcu*fczu,co);
%      ph=plot3(xcu,ycu,zcu*fczu,col(1,jj));
      set(ph,'linewidth',ilin)
      if shav==6 | shav==8
       set(ph,'linewidth',1)
      end
%      end
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
   end
%   ' map' , keyboard
