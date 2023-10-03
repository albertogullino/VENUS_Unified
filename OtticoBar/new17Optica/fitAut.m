function ve=vefit(Vei,x,x0,gr)
Ve=log10(abs(Vei));
s=size(Ve,1);
iver=0;
if iver==1
 h=figure;
end
for k=1:s
 coA=polyfit(x,Ve(k,:),gr);
 ve(k,1)=polyval(coA,x0);
 vev=polyval(coA,x);
 if iver==1
  figure(h)
  semilogy(x,abs(Ve(k,:)),x0,abs(ve(k)),'wo',x,abs(vev),'r.'), pausak
 end 
end

ve=10.^ve;
%figure, semilogy(abs(Ve)), hold on, semilogy(abs(ve),'w.'), pausak