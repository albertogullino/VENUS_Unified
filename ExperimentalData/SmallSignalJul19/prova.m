f=logspace(-1,1,101);
o=2*pi*f*1e9;
jo=1i*o;

Z0=50;
Cp=20e-15;

Rs=7;
Ra=90;
Ca=600e-15;

Ga=1/Ra;

Y=jo*Cp+(Ga+jo*Ca)./(1+Rs*(Ga+jo*Ca));

figure, plot(f,real(Y),f,imag(Y)), pausak

z=1./(Y*Z0);
S=(1-Y*Z0)./(1+Y*Z0);

figure
        P=polar(angle(S),abs(S));
        set(P,'Linewidth',2)
        pausak
        
            figure, plot(f,real(z),f,imag(z))