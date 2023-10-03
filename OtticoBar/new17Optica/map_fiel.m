function map_fiel(XP,YP,Emoto,aax,xcu,ycu,zcu,titl,xcum,ycum,ibar)

    surf(XP,YP,Emoto),
    shading('interp'),
    view(0,90),
    axis square,
    axis equal,
    axis([-1 1 -1 1]*aax),
    axis off
    fczu=1;
    title(titl),
    hold on, fill3(xcu,ycu,zcu*fczu,'w')
    if length(xcum)==length(zcu)
     fill3(xcum,ycum,zcu*fczu,'y')
    end
    if ibar==1
     colorbar
    end
