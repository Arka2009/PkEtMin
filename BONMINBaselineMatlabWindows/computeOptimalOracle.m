function optx = computeOptimalOracle(D,benchid,AET,BET,AP,BP,LLIM,M,N)
    % Compute the oracle, with given deadline, benchmark composition 
    % benchmark characteristics (AET,BET,AP,BP,LLIM) and number of 
    % phases (M) and maximum number of cores in the system (N)

    % Set the objective function
    objfunc = @(x) (x(M+1));

    % Setup the constraints
    pkpFn = @(a,b,pkp) (computePowerPerPhase(a,b,AP,BP,LLIM,M,N) - pkp);
    constraintsfunc = @(x) [arrayfun(pkpFn,x(1:M),transpose(benchid),x(M+1)*ones(M,1));computeExecTime(x(1:M),transpose(benchid),AET,BET,LLIM,M,N)];
    cl = [transpose(-inf*ones(1,M));0];
    cu = [transpose(zeros(1,M));D];

    % Variable Bounds
    lb = [transpose(LLIM(benchid));0];
    ub = [transpose(N*ones(1,M));inf];

    % Integer Constraints
    xtype = join([char('I'*ones(1,M)),'C']);
    x0    = [N*ones(M,1);0.0];

    % Create OPTI Object
    % bonminopts = bonminset('algorithm','B-Hyb');
    maxtime = 1.05*(10^(M/10));
    opts = optiset('solver','bonmin','display','iter',...
    'maxiter',1e8,'maxtime',maxtime);

    Opt  = opti('fun',objfunc,'nl',constraintsfunc,cl,cu,'bounds',lb,ub,...
    'x0',x0,'options',opts,'xtype',xtype);
    

    % Solve the MINLP problem
    [x,fval,exitflag,info] = solve(Opt,x0);
    optx = x;
    % optx = N*ones(M,1);
end