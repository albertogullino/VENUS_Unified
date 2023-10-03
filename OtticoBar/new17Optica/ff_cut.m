global fafir xc1f
Ef2=Ef.^2;

liman=89;
fi0=find(abs((fian0-pi))<1e-4);
xc1=flipud(X(2:end,fi0));
ff1=flipud(Ef2(2:end,fi0));
maff1=Ef2(1,fi0);
fi0=find(fian0==0);
xc1=[xc1; X(2:end,fi0)];
ff1=[ff1; Ef2(2:end,fi0)];
xc1f=linspace(-liman,liman,200);
ff1f=spline(xc1,ff1/maff1,xc1f);
%figure, plot(xc1,ff1/maff1,'.',xc1f,ff1f), pausak
fi0=find(abs((fian0-pi/2*3))<1e-4);
xc2=flipud(Y(2:end,fi0));
ff2=flipud(Ef2(2:end,fi0));
maff2=Ef2(1,fi0);
fi0=find(abs((fian0-pi/2))<1e-4);
xc2=[xc2; Y(2:end,fi0)];
ff2=[ff2; Ef2(2:end,fi0)];
xc2f=linspace(-liman,liman,200);
ff2f=spline(xc2,ff2/maff2,xc2f);
%figure, plot(xc1,ff1/maff1,'g*',xc2,ff2/maff2,'r*'),
%fafir=ff1f/max(ff1f);
fafir(:,1)=ff1f'/max(ff1f);
fafir(:,2)=ff2f'/max(ff2f);

 if ifp==-10
  figure, plot(xc1f,ff1f/max(ff1f),'r',xc2f,ff2f/max(ff2f),'g'),
  axis([-liman liman 0 1])
  pausak
 end

%Please find attached far-field measurements on device 5_90 at 1mA in
%[01-1]-direction (measurement A) and [011]-direction (measurement B).
%Column one is far-field angle in degrees, column two the normalized
%far-field intensity. I hope they agree with your simulations; as you can
%see, the main difference is the appearance of diffraction side lobes
%when I measure along the minor ellipse axis. I also have measurements at
%higher currents but in the multi-mode range, the curves get somewhat
%messy. Let me know if you need anything else.
%
%> P.S. Could you please tell me where the short manuscript on the relief structure  was sent for publication and if it is published yet?
%
%It was published in Electronics Letters, vol.38, no.2, pp.77-78, 2002


