load sav

fir=2
nr1=nv(fir,1);
nr2=nv(fir,2);
%ifpT=-10;

Kr=KK(15);

ifpT=5
mb=[1 3 5];             
mb=[1 3 5];             
%mb=[1 ];             
%mb=[1 3 5]-1;             
%             [Oo1,Oo2,G1c,G2c]=Teq1(Kr,lar,dos,period,DC,nr1,nr2,r_in,r_out,rr,mb,ifpT);       

Kr=0

             [Oo1c,Oo2c,G1c,G2c]=Teq1p(Kr,lar,dos,period,DC,nr1,nr2,r_in,r_out,rr,mb,ifpT);       

map(abs(Oo1c)), shading faceted          , pausak
map(abs(Oo2c)), shading faceted          
%             CONT1=sum(sum(abs(Oo1-Oo1c)))
%             CONT2=sum(sum(abs(Oo2-Oo2c)))