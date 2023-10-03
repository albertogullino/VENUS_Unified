Tp=1;
for kma=40:56
 Tp=Tv(:,:,kma)*Tp;
end
Tup=Tp;

Tp=1;
for kma=4:5
 Tp=Tv(:,:,kma)*Tp;
end
Tp1=Tp;

Tp=1;
for kma=6:22
 Tp=Tv(:,:,kma)*Tp;
end
Tp20=Tp;
Tp2=Tp^17;

Tp=1;
for kma=23:38
 Tp=Tv(:,:,kma)*Tp;
end
Tp3=Tp;

Tmir=Tp3*Tp2*Tp1;