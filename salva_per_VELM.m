if indVELM>=1
    %salvo input di VELM
    VELMInput(indVELM).matgain=mode.matgain;
    VELMInput(indVELM).DeltaN=mode.DeltaN;
    VELMInput(indVELM).Deltan=mode.Deltan;
    VELMInput(indVELM).efield_y=mode.efield_y;
    VELMInput(indVELM).efield_x=mode.efield_x;
    VELMInput(indVELM).elecABS=mode.elecABSvelm;
    VELMInput(indVELM).holeABS=mode.holeABSvelm;
    VELMInput(indVELM).vlambda=mode.vlambda;
    VELMInput(indVELM).alpha=mode.alpha;
    VELMInput(indVELM).Lm=mode.Lm;
    VELMInput(indVELM).NQW=mode.NQW;
    VELMInput(indVELM).Gamma_z=mode.Gamma_z;
    VELMInput(indVELM).nindexQW=mode.nindexQW;
    VELMInput(indVELM).fPES=mode.fPES;
    VELMInput(indVELM).fPdif=mode.fPdif;
    VELMInput(indVELM).E2=mode.E2;
    
    VELMInput(indVELM).TmVelm=mode1.TmVelm;
    VELMInput(indVELM).LamVelm=mode1.LamVelm;
    VELMInput(indVELM).DeltaT=mesh.DeltaT;
    VELMInput(indVELM).DeltaTvelm=mesh.DeltaTvelm;
    VELMInput(indVELM).ygrid=mesh.ygrid;
    VELMInput(indVELM).xgrid=mesh.xgrid;
    VELMInput(indVELM).nnx=mesh.nnx;
    VELMInput(indVELM).nny=mesh.nny;
    VELMInput(indVELM).v0_dd=v0_dd;
    VELMInput(indVELM).DeltaTmax=DeltaTmax;
    VELMInput(indVELM).indVoltage=indv;
    VELMInput(indVELM).indVELM=indVELM;
    VELMInput(indVELM).uvet=uvet;
    
    if indVELM==1
        VELMInput(indVELM).Gmod=[];
    else
        VELMInput(indVELM).Gmod=mode.Gmod;
        
        % save in VelmSa.mat
        eval(['save ',VELMinput,' VELMInfo VELMInput ParVet VelmOptions mode mode1 mesh fil_str vind'])
    end
end