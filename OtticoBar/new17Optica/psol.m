imi=find(fap<=0)+1;
%imi
%keyboard
if length(imi)>0
imi=imi(1);
if imi==length(fap)
 imi=imi-1;
end
if imi>1  & imi<nl

fi=imi+[-1 0 1];
cof=polyfit(dlav(fi),fa_tot(fi),2);
dlavd=roots(cof);
[du,iro]=min(abs(dlavd-dlav(imi)));
dlav0=dlavd(iro);
cog=polyfit(dlav(fi),log10(gth(fi)),2);
gth0=10^polyval(cog,dlav0);
if isnan(gth0)==1
cog=polyfit(dlav(fi),(gth(fi)),2);
gth0=polyval(cog,dlav0);
end

coga=polyfit(dlav(fi),log10(fve(fi)),2);
gam0=10^polyval(coga,dlav0);
coen=polyfit(dlav(fi),log10(eneu(fi)),2);
en0=10^polyval(coen,dlav0);
if gth0>0
 iso=iso+1;
 enmv(iso)=real(en0); 
 dlv(iso)=real(dlav0);
 gtv(iso)=real(gth0);
 gamv(iso)=real(gam0); 
 dlv0(iso)=(isol-1)*2*pi;
% 'isol ', isol, pausak
end
end %if
end %length