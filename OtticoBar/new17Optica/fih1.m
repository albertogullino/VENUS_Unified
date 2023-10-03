load saca

figure, plot(zet,real(Ez),'r.',zet,imag(Ez),'c.'), pausak
figure, plot(zet,abs(Ez).^2,'r.'), pausak

I=abs(Ez).^2;
ia=20;
ib=21;
fa=I(ia);
fap=diff(I(ia+[-1 0]))/diff(zet(ia+[-1 0]));

fb=I(ib);
fbp=diff(I(ib+[0 1]))/diff(zet(ib+[0 1]));
d=diff(zet([ia ib]));

za=zet(ia);
zb=zet(ib);

zeh=linspace(za,zb,50);

df=fb-fa;
k=fbp/df;

fap=fa+df*exp(k*(zeh-zb));

figure, plot(zet,I,'r.',zeh,fap,'g.'),
return
