load sav

fir=2
nr1=nv(fir,1);
nr2=nv(fir,2);
%ifpT=-10;

Kr=KK([1]);

mb=1;
             [Oo1,Oo2,G1,G2]=Teq(Kr,lar,dos,period,DC,nr1,nr2,r_in,r_out,rr,mb,ifpT);   

             [Oo1c,Oo2c,G1c,G2c]=Teq1(Kr,lar,dos,period,DC,nr1,nr2,r_in,r_out,rr,mb,ifpT);       
             
             CONT1=sum(sum(abs(Oo1-Oo1c)))
             CONT2=sum(sum(abs(Oo2-Oo2c)))