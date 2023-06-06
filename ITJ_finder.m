iTJ=[];

for iNTJ=1:mode.nBTJ % = size(mesh.IBTJ,1)
    for itj=1:size(mesh.IBTJ,2)
        indTJ=mesh.IBTJ{iNTJ,itj}(1:length(mesh.IBTJ{iNTJ,itj})/mesh.nnx);
        iTJ=[iTJ, indTJ];
    end
    iTJ=unique(sort(iTJ));
end

% iTJ;

IBTJ=unique(sort(cell2mat(mesh.IBTJ)));

% meshPLOT
