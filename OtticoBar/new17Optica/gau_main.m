clear
close all
a=5;
wv=logspace(0,2,11);

for kw=1:length(wv)
 w=wv(kw)
 gau_sub
 Re(kw)=R_eq;
 'fine apertura', keyboard
end

figure, semilogx(wv,Re)
pausak
figure, semilogx(wv,-log10(1-Re))