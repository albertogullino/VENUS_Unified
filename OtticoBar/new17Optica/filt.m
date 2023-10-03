function f=filt(y,so)

me=mean(abs(y));
ma=max(abs(y));
f=y;
if ma/me>.3
 fin=find(y>so*me);
 f(fin)=0;
end
