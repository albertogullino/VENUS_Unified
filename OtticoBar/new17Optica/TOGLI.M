ayis=ayi;
da=diff(ayis);
so=.05;
fiv=find(da>so);
fivn=find(da<=so);
ayu=ayi;
ayu(fivn)=NaN;
figure, plot(ayi,'.'), hold on, plot(ayu,'ro')
