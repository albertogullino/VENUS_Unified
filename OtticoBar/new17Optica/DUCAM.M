
    E2xp=(Exdu*ensx);
    E2yp=(Eydu*ensy);
    E2zp=(Ezdu*ensz)*zete;

    figure,
       set(gcf,'Position',pograd);

       aax=5;

       ibar=1;
       iaoff=0

       titl=' E_x  Output';
       subplot(1,3,1)
        map_fnew(XP,YP,E2xp,aax,Cug.x,Cug.y,Cug.z,titl,ibar,iaoff)
       titl=' E_y  Output';
       subplot(1,3,2)
        map_fnew(XP,YP,E2yp,aax,Cug.x,Cug.y,Cug.z,titl,ibar,iaoff)
       titl=' E_z  Output';
       subplot(1,3,3)
        map_fnew(XP,YP,E2zp,aax,Cug.x,Cug.y,Cug.z,titl,ibar,iaoff)
        pausak
