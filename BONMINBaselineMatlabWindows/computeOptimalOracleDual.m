function optx = computeOptimalOracleDual(P,benchid,AET,BET,AP,BP,LLIM,M,N)
    % Compute the oracle, with given peak power, benchmark composition 
    % benchmark characteristics (AET,BET,AP,BP,LLIM) and number of 
    % phases (M) and maximum number of cores in the system (N)

    % Set the objective function
    objfunc = @(x) (computeExecTime(x,transpose(benchid),AET,BET,LLIM,M,N));

    % Setup the constraints
    pkpFn = @(a,b) (computePowerPerPhase(a,b,AP,BP,LLIM,M,N));
    constraintsfunc = @(x) [arrayfun(pkpFn,x,transpose(benchid))];
    cl = zeros(M,1);
    cu = P*ones(M,1);

    % Variable Bounds
    lb = transpose(LLIM(benchid));
    ub = transpose(N*ones(1,M));

    % Integer Constraints
    xtype = char('I'*ones(1,M));
    x0    = N*ones(M,1);

    % Create OPTI Object
    opts = optiset('solver','bonmin','display','iter');
    Opt = opti('fun',objfunc,'nl',constraintsfunc,cl,cu,'bounds',lb,ub,...
    'x0',x0,'options',opts,'xtype',xtype);

    % Solve the MINLP problem
    [x,fval,exitflag,info] = solve(Opt,x0);
    optx = x;
end