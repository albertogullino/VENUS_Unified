% PnQW=
function [nQW,pQW]=Extract_carrierQW(mode)

for iV=1:length(mode.ii_dd)
%     for iQW=1:P.mesh.nnxQW{1}
        Niv=mode.nQW(iV);
        nQW(iV,:)=Niv{1}{2};
        
        Piv=mode.pQW(iV);
        pQW(iV,:)=Piv{1}{2};
end