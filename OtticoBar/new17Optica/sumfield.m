%'sf 10', keyboard
if iant==0
 Andup=Andu.*pes(Pus);
 AnduHp=AnduH.*pes(Pus);
else
 Andup=Andu;
 AnduHp=AnduH; 
end

%'quin ando', keyboard
  if iLP==0

   AnFs=Andup.*Fs;

%   Mvefm=diag(AnFs)*[ Mvefm0;  segem*Mvefm0];
%   Mvegm=diag(AnFs)*[ Mvegm0;  segem*Mvegm0];
%
%   Mvefp=diag(AnFs)*[ Mvefp0; -segem*Mvefp0];
%   Mvegp=diag(AnFs)*[-Mvegp0;  segem*Mvegp0];

   Mvefm=diag(AnFs)* Mvefm0;
   Mvegm=diag(AnFs)*Mvegm0;

   Mvefp=diag(AnFs)*Mvefp0;
   Mvegp=diag(AnFs)*Mvegp0;
   
   if iztm==1
    AnFsz=j*2*fatEz*AnFs.*fkt.*Fs;
    Mvez=diag(AnFsz)*Mvez0;
   end  



   
   Exdu=besp'*Mvefp+besm'*Mvefm;
%   Eydu=-besp'*Mvegp+besm'*Mvegm;
%   'corretto '
   Eydu=besp'*Mvegp+besm'*Mvegm;

   AnFs=AnduHp.*Fs;

   Mvefm=diag(AnFs)* Mvefm0;
   Mvegm=diag(AnFs)*Mvegm0;
   Mvefp=diag(AnFs)*Mvefp0;
   Mvegp=diag(AnFs)*Mvegp0;
   
   Hydu=besp'*Mvefp+besm'*Mvefm;
   Hxdu=-(besp'*Mvegp+besm'*Mvegm);

   if iztm==1
    Ezdu=besz'*Mvez;
   else
    Ezdu=0;
   end

%   'corretto contr ', keyboard
  else

   AnFs=Andup.*Fs;
   Mvefm=diag(AnFs)*Mvefm0;
   Exdu=besm'*Mvefm;
   Eydu=0;
   Ezdu=0;
   AnFsH=AnduHp.*Fs;
   MvefmH=diag(AnFsH)*Mvefm0;
   Hydu=besm'*MvefmH;
   Hxdu=0;

  end

%     if fimaxi<2*pi
%      Exdu=[Exdu fliplr(Exdu) ];
%      if iLP==0
%       Eydu=[Eydu fliplr(Eydu) ];
%       Ezdu=[Ezdu fliplr(Ezdu) ];
%      end
%     end

     if fimaxi<2*pi
      lfil=length(fian)-1;
      Exdu=[Exdu(:,1:lfil) fliplr(Exdu) ];
      if iLP==0
       Eydu=[Eydu(:,1:lfil) fliplr(Eydu) ];
       Ezdu=[Ezdu(:,1:lfil) fliplr(Ezdu) ];
      end
     end

%     if fimaxi<2*pi
%      Exdu=[Exdu fliplr(Exdu) Exdu(:,1)];
%      if iLP==0
%       Eydu=[Eydu fliplr(Eydu) Eydu(:,1)];
%       Ezdu=[Ezdu fliplr(Ezdu) Ezdu(:,1)];
%      end
%%     end
%' passo sumfield ', keyboard
