  if iLP==0

   AnFs=AnduH.*Fs./Znor(Pus);
   AnFsz=-2*j*Andu.*fkth.*Fs./Znor(Pus);

   Mvefm=diag(AnFs)*[Mvefm0; Mvefm0];
   Mvegm=diag(AnFs)*[Mvegm0; Mvegm0];

   Mvefp=diag(AnFs)*[Mvefp0; Mvefp0];
   Mvegp=diag(AnFs)*[Mvegp0; Mvegp0];

   if iztm==1
    Mvez=diag(AnFsz)*[Mvhz0; Mvhz0];
   end

   Hydu=besp'*Mvefp+besm'*Mvefm;
   Hxdu=-(-besp'*Mvegp+besm'*Mvegm);


   if iztm==1
    Hzdu=besz'*Mvez;
   else
    Hzdu=0;
   end

  else

   AnFs=Andu.*Fs;
   Mvefm=diag(AnFs)*Mvefm0;
   Hydu=besm'*Mvefm;
   Hxdu=0;
   Hzdu=0;

  end


     if fimaxi<2*pi
      lfil=length(fian)-1;
      Hxdu=[Hxdu(:,1:lfil) fliplr(Hxdu) ];
      if iLP==0
       Hydu=[Hydu(:,1:lfil) fliplr(Hydu) ];
       Hzdu=[Hzdu(:,1:lfil) fliplr(Hzdu) ];
      end
     end

