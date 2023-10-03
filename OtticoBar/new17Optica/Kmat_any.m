%' kmat_any', keyboard
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%/%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Contributo dell'Anisotropia  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if pasnu==2

 if isfield(Pf,'iany_CIRC')==1
  iany_CIRC=Pf.iany_CIRC;
 else
  iany_CIRC=0;
 end

 kmat_any_pasnu2New
 
else

 %kmat_any_pasnu1
 kmat_any_pasnu3

end


%'fine aby', keyboard
 