function Sint=spacint(F,rdr,lfp,pesf,Nsett)

Fi=F(:,1:lfp);
Sint=pesf*(Fi'*rdr)*Nsett;
%disp(' spacint '), keyboard
