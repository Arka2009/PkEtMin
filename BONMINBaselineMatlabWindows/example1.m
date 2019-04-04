% Objective
fun = @(x) -x(1) - x(2) - x(3);    

% Linear Constraints
A = [1 -1 0 0;
     1  0 1 1];
b = [0;2];

% Nonlinear Constraint
nlcon = @(x) (x(2) - 0.5)^2 + (x(3) - 0.5)^2;
nlrhs = 0.25;
nle = -1; % -1 for <=, 0 for ==, +1 >=        

% Bounds
lb = [0;0;0;0];
ub = [1;Inf;Inf;5];

% Integer Constraints
xtype = 'BCCI';

%Initial Guess
x0 = [0;0;0;0];

% Create OPTI Object
Opt = opti('fun',fun,'nlmix',nlcon,nlrhs,nle,'ineq',A,b,'bounds',lb,ub,...
           'xtype',xtype)

% Solve the MINLP problem
[x,fval,exitflag,info] = solve(Opt,x0)