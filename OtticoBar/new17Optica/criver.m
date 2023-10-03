load sa

 Ga1_old=Ga1;

% [Gacrit,Trcrit,Trc,Grc]=Gam_critUn(Tstor,Ga1,Mcrit,ficri,fmlsi);
 [Gacrit,Trcrit,Trc,Grc]=Gam_critScatt(Ga1,Mcrit,ficri,fmlsi);
% ' qui cont nuovo dopo', keyboard

% metodo solito
Tet=1;
for kf=ficri'
 Ted=expm(Mcrit{kf});
 Te{kf}=Ted;
 Tet=Ted*Tet;
end

p1=1:length(KKt);
p2=p1+length(KKt);

T11=Tet(p1,p1);
T12=Tet(p1,p2);
T21=Tet(p2,p1);
T22=Tet(p2,p2);
Geq=(T11*Ga1_old+T12)*inv(T21*Ga1_old+T22);

map(log10(abs(Geq)))
pausak
map(log10(abs(Gacrit)))
pausak
map(abs(Gacrit-Geq))
