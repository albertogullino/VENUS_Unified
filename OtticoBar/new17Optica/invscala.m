      Anout=[Afzf(:,iFFsect); Afz(:,iFFsect)];
    %  Afzf_s=Afzf;
    %  Afz_s=Afzc;
%save gaussOo OoScala
load gaussOo
Oo=inv(OoScala);
      Anout_prima=Oo*Anout;
      Afzf(:,iFFsect)=Anouted(1:end/2);
      Afz(:,iFFsect)=Anouted(end/2+1:end);
      Anoute=Anouted(1:end/2)+Anouted(end/2+1:end);
      Anout=Anoute;
