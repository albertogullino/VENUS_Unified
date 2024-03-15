function Line=line_shape(delE,Gamma0,ganm,iline)

% All quantities are returned in 1/s

switch iline
    
    case 'lorentzian'
        % Classic lorentzian profile; see for example 1988AhnD_JQE, eq. (8)
        Line=delE-1i*Gamma0;
    case 'nMark'
        % Lineshape including non-Markovian effects, taken from
        % 1991Tomita_JQE, "A new density matrix theory for semiconductor
        % lasers, including non-Markovian intraband relaxation and its
        % application to nonlinear gain", eq. (44)
        %$
        A=(delE)/ganm;
        B=-1i*Gamma0/ganm;
        N=max([5 fix(abs(max(max(B))))])*3;
        N=min([N,20]);
        sa=A+B;
        l=1./(sa+N);
        for in=N-1:-1:0
            l=(1+B.*l)./(sa+in);
        end
        Line=ganm./l;
    otherwise
        error('Error in line_shape.m ')
end
