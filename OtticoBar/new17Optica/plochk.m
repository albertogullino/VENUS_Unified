load plo

fim10=find(Fpl==-10);
fim5=find(Fpl==-5);
fi0=find(Fpl==0);
fim20=find(Fpl==-20);
fima=find(Fpl>0);
td=Tpl(fim10);
fd=Fpl(fim10);
T=td;
F=fd;
ltd=length(td);
if ltd>2
 T=reshape(td,2,ltd/2);
 F=reshape(fd,2,ltd/2);
 T=[T; T(1,:)*NaN];
 F=[F; T(1,:)*NaN];
 ltd=prod(size(T));
 T=reshape(T,ltd,1);
 F=reshape(F,ltd,1);
end
T10=T;
F10=F;

td=Tpl(fim5);
fd=Fpl(fim5);
T=td;
F=fd;
ltd=length(td);
if ltd>2
 T=reshape(td,2,ltd/2);
 F=reshape(fd,2,ltd/2);
 T=[T; T(1,:)*NaN];
 F=[F; T(1,:)*NaN];
 ltd=prod(size(T));
 T=reshape(T,ltd,1);
 F=reshape(F,ltd,1);
end
T5=T;
F5=F;

td=Tpl(fi0);
fd=Fpl(fi0);
T=td;
F=fd;
ltd=length(td);
if ltd>2
 T=reshape(td,2,ltd/2);
 F=reshape(fd,2,ltd/2);
 T=[T; T(1,:)*NaN];
 F=[F; T(1,:)*NaN];
 ltd=prod(size(T));
 T=reshape(T,ltd,1);
 F=reshape(F,ltd,1);
end
T0=T;
F0=F;

td=Tpl(fim20);
fd=Fpl(fim20);
T=td;
F=fd;
ltd=length(td);
if ltd>2
 T=reshape(td,2,ltd/2);
 F=reshape(fd,2,ltd/2);
 T=[T; T(1,:)*NaN];
 F=[F; T(1,:)*NaN];
 ltd=prod(size(T));
 T=reshape(T,ltd,1);
 F=reshape(F,ltd,1);
end
T20=T;
F20=F;
if max(Tpl)<500
 tat=1;
 xl='elapsed time (s)';
elseif max(Tpl)>=500  &  max(Tpl)<10000
 tat=1/60;
 xl='elapsed time (min)';
else
 tat=1/3600;
 xl='elapsed time (hour)';
end

iskif=0;
if iskif==0
   figure
   plot(tat*T10,F10,'r.-',tat*T5,F5,'c.-',...
   tat*T20,F20,'b.-',tat*T0,F0,'g.-',tat*Tpl(fima),Fpl(fima),'w.-')

axis([0 max(Tpl)*tat -25 max(Fpl)]), grid,
legend(' -10 Coup. Matrix',' -5 prod sopra ',' -20 Eig ',' 0 prod sotto ',2)
xlabel(xl)
end
