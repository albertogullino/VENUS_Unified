load Janoptics
ES=squeeze(Es)/max(Es);
figure, plot(lambdavet*1e9+1.8,ES,'linewidth',2)
hold on
load PL
plot(PL(:,1),PL(:,2),'o')
ylabel('PL, arb.un.')
xlabel('Wavelength, nm')
xlim([780 880])