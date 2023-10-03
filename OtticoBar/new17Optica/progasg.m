  T1=Tdu;
  T2=Tab_air;
  Ta=T_air;
  T=T1*T_air*Tab_air;  
  dGas=diag(Gas);
  s=size(Tdu);
  l1=s(1)/2;
  l2=s(1)/2+1;
  l3=s(1);
  T11=T(1:l1,1:l1);
  T12=T(1:l1,l2:l3);
  T21=T(l2:l3,1:l1);
  T22=T(l2:l3,l2:l3);
  Geq=(T11*dGas+T12)*inv(T21*dGas+T22);
  
  
    T11=T2(1:l1,1:l1);
    T12=T2(1:l1,l2:l3);
    T21=T2(l2:l3,1:l1);
    T22=T2(l2:l3,l2:l3);
  Geq1i=(T11*dGas+T12)*inv(T21*dGas+T22);
  
  T1a=T1*Ta;
  
      T11=T1a(1:l1,1:l1);
      T12=T1a(1:l1,l2:l3);
      T21=T1a(l2:l3,1:l1);
      T22=T1a(l2:l3,l2:l3);
  Geq1=(T11*Geq1i+T12)*inv(T21*Geq1i+T22);
  
map(Geq1-Geq)
  
  