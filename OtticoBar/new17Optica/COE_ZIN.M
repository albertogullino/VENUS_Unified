izok=0;
fical=1;
if length(Lizi)>1
 Lzdu=Lizi-zi(end);
 fima=find(Lzdu>0);
% [dos Lzdu]
% pausak
 if length(fima)>0 & min(Lzdu(fima))<dos
  izok=1;
  fical=find(Lzdu<dos & Lzdu>0);
  Lizv=sort([0 Lzdu(fical)]);
  Liz=diff(Lizv([0 1]+1));
 else
  Liz=100;
 end
else
  Liz=100;
end

%if dos>0 & izok==1
if dos>0
ioxi=0;
ncompari=0;
     Oos=Oo;
     Dos=dos;
%     'inizio coe_zi: Dos=, istr=',
%     [Dos istr],
     istr
%     Dos
%     if istr==32
%      'istr=32'
%      keyboard
%     end
%     pausak
     Amods=Amod;
     zip=zi(icfiez);

     if Dos>Liz & ifiez==1
%      ncompari=nitr(1)+10*shtr(1)
      ncompari=nitr(1)+10*aitr(1);
%      if ncompari>10
%       'ossido', keyboard
%       ioxi=1;
%      end
      finc=find(abs(ncompari-ncompar)<1e-6);
%       'valuto Oo ', keyboard
      if length(finc)==0

%       pausak
       ncompar(incdiv)=ncompari;
       
       eval(emme)
       Treus{incdiv}=Oo;
       incdiv=incdiv+1;
      else
       Oo=Treus{finc};
      end
      icfiez=icfiez+1;
      zi(icfiez)=zi(icfiez-1)+Liz;
%      ' uno ', keyboard
      Nz(icfiez)=nitr(1);
%      ' uno '
%      pausak
      Amod=Oo*Amod;
      Acoz(:,icfiez)=Amod;
      Lis=2*Liz;
      ifirst=0;
      izemp=0;
%      if Dos>10
%       ' DOs'
%       Dos=100;
%%       keyboard
%      end
%      tic
%  Lis
%  Dos
% 'Lis ', pausak
%      while Lis<Dos

      decis=Lis-Dos;
%      keyboard

%      while decis<0
     for kkz=2:length(fical)
%       ' dentro while Lis', pausak
       icfiez=icfiez+1;
%       'altro'
%       'sue ', keyboard
       Liz=diff(Lizv([0 1]+kkz));
       zi(icfiez)=zi(icfiez-1)+Liz;
       Nz(icfiez)=nitr(1);
%       pausak
%%       if izemp==0
        Amod=Oo*Amod;
%%       else
%%        Amod(pu1)=O11.*Amod(pu1)+O12.*Amod(pu2);
%%        Amod(pu2)=O21.*Amod(pu1)+O22.*Amod(pu2);
%        ' proviamo ', keyboard
%        'vel'
%%       end
       Acoz(:,icfiez)=Amod;
       Lis=Lis+Liz;
       decis= Lis-Dos;
       if abs(Lis-Dos)<=1e-6 | decis>0
        decis=0;
       end
%       ' decis '
%       decis
%       pausak
      end
%      toc
%      keyboard
     end
     icfiez=icfiez+1;
%     'altro 1'
%     'tre', keyboard
      Nz(icfiez)=nitr(1);
     zi(icfiez)=zip+Dos;
%     pausak
     Amod=Oos*Amods;
     Acoz(:,icfiez)=Amod;

%      zi(end-2:end)
%     'fine coe_zi', pausak
%if dos<1e-5
% ' controlla coezi ', keyboard
%end
end
%if istr==48
% ' controlla coezi NNN ', keyboard
%end
