%'entro mod_ff', keyboard
if iFFsect==1
 rff=real(nv(1,1));
else
 rff=real(nv(end,1));
end
if rff<1
 rff=1;
end
rrfa=real(rr/rff);
 if length(lfi_inp)<2
  nrff=255;
 else
  nrff=lfi_inp(2);
 end
 

 iR=find(KK*rrfa<1);
 iRm=iR(end);
 Kff=linspace(.0,KK(iRm),nrff)';
 if iRm<length(KK)
  Kff=linspace(0.0,1/rrfa,nrff)';
 end

 %nrff=355;
 %Kff=linspace(0,.1,nrff)';

% iR=find(Kff*rr<=1);
% iRm=iR(end);
% teRma=asin(Kff(iRm)*rr);
 teRma=asin(Kff(end)*rrfa);
% ' Kff ', keyboard

% iR=find(KK*rr<1);
% iRm=iR(end);
% teRma=asin(KK(iRm)*rr);
% teR=linspace(0,teRma,nrff)';
 teR=asin(Kff*rrfa);
 X1=teR*180/pi*cos(fian);
 Y1=teR*180/pi*sin(fian);
 Xpr=sin([teR])*cos(fian);
 Ypr=sin([teR])*sin(fian);
 Zpr=cos([teR])*ones(size(fian));

%'% calcolo Im ', keyboard
intgu_ff

npkf=nrff;
