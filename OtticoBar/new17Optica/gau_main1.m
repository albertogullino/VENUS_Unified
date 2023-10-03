clear
close all
a=5;
wv=1*logspace(0,2,11);

for kw=1:length(wv)
 w=wv(kw)
 w2=w^2;
 gauB_su
 Re1(kw)=R1_eq;
 Re2(kw)=R2_eq;
 Te1(kw)=T1_eq;
 Te2(kw)=T2_eq; 
 'fine apertura cilindrica', keyboard 
end

figure, semilogx(wv,-log10(1-Re1),wv,-log10(1-Re2)), pausak
figure, semilogx(wv,Re1,wv,Re2)
figure, loglog(wv,Te1,wv,Te2)