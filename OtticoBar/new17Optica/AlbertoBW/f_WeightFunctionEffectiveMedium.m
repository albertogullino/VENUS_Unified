
function Objective=f_WeightFunctionEffectiveMedium(n2,S11,n1,n3,L,lambda,theta,Polarization)

[S11test,S21test] = f_ComputeSlab(n2,n1,n3,L,lambda,theta,Polarization);

Objective=log10(abs(S11test-S11));
% Objective=abs(S11test)-abs(S11); % non è l'ideale, ripensarci
% Objective=real(S11test)-real(S11); % non è di nuovo l'ideale, ripensarci

% figure(9432)
% hold on
% plot(n2,log10(abs(Objective)),'*'),drawnow

return