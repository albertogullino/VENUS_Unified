function [nQW,pQW,Efn2D,Efp2D,Ccapn,Ccapp]=Extract_carrierQWcell(mode)

for iV=1:length(mode.ii_dd)
        % Bound carriers
        Niv=mode.nQW(iV);   % mode.nQW: vector of cells, long as the bias vector
        nQW(iV,:)=Niv{1}{2};    % first index access to the cell, the second is the QW index!
        
        Piv=mode.pQW(iV);
        pQW(iV,:)=Piv{1}{2};
        
%         % Quasi-Fermi levels
%         Efiv=mode.EFn2D(2,iV); % matrix: [indQW,indBias]
%         Efn2D(iV,:)=Efiv{1};    % the index is needed to access to the cell
%         
%         Efiv=mode.EFp2D(2,iV);
%         Efp2D(iV,:)=Efiv{1};
        
end

% Quasi-Fermi levels
Efn2D=cell2mat(mode.EFn2D(2,:)');   % the index is needed to access to the cell
Efp2D=cell2mat(mode.EFp2D(2,:)');   % the index is needed to access to the cell


% Ccap (mode.Ccapn has the same dimension of mode.Efn2D!!!)
Ccapn=cell2mat(mode.Ccapn(2,:)');   % the index is needed to access to the cell
Ccapp=cell2mat(mode.Ccapp(2,:)');   % the index is needed to access to the cell


