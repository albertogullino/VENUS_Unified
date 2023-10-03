%po_pru
%if ipolar~=0
if pasnu~=1
% po_new
 po_new_all
 
else
 po_new_all
end

 Xpr=sin([teR])*cos(fian);
 Ypr=sin([teR])*sin(fian);
 Zpr=cos([teR])*ones(size(fian));
 Pz_n=Pvectz.*Zpr;
% if ifp==-10
% ' controllo Point ', keyboard
%   cont_ff
%   pausak
% end

 if iLP==0
 Px_n=Pvectx.*Xpr;
 Py_n=Pvecty.*Ypr;
 Poim_old=sqrt(Pvectx.^2+Pvecty.^2+Pvectz.^2);
 Poim_n=sqrt(Px_n.^2+Py_n.^2+Pz_n.^2);
  if iFFte==1
   Poim=Poim_n;
  else
  
   Poim=Poim_old;
  end
 else
 Poim=abs(Pvectz);
 end
 fi0=find(Poim==0);
 Poim(fi0)=1e-50;
 Pointing=Poim;
 te=acos(Pvectz./Poim);

 if iFFte==1
  X=X1;
  Y=Y1;
 else
  X=rox*fm0;
  Y=rox*gm0;
 end
 Ef=sqrt(Pointing);

