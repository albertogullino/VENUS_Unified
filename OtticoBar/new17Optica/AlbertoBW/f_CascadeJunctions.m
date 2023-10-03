
function [S11,S12,S21,S22]=f_CascadeJunctions(JunctionInfo,indJunc1,indJunc2)

for indJunction=indJunc1:indJunc2
    if indJunction==indJunc1
        S11p=JunctionInfo(indJunction).S11;
        S21p=JunctionInfo(indJunction).S21;
        S12p=JunctionInfo(indJunction).S12;
        S22p=JunctionInfo(indJunction).S22;
    else
        S11s=JunctionInfo(indJunction).S11;
        S21s=JunctionInfo(indJunction).S21;
        S12s=JunctionInfo(indJunction).S12;
        S22s=JunctionInfo(indJunction).S22;
        kzL=JunctionInfo(indJunction).kzL;
        deltaz=JunctionInfo(indJunction).zj-JunctionInfo(indJunction-1).zj;
        [S11p,S12p,S21p,S22p]=f_CascS(S11p,S12p,S21p,S22p,S11s,S12s,S21s,S22s,kzL,deltaz,inf);
    end
end

S11=S11p;
S21=S21p;
S12=S12p;
S22=S22p;

return