function [xstar,fxstar,k,exitflag,xsequence] = myfmincon(fun,x0,A,b,C,d,p,q,myoptions)
% MYFMINCON Attempts to solve the problem:
%                   min f(x)  
%                     s.t.
%                   A*x = b
%                  C*x >= d
%                  g(x) = 0
%                 h(x) >= 0
% and, if successful, returns a local minimizer xstar and the related local
% optimum fxstar=f(xstar). The solver employs a Sequential Quadratic programming
% optimization scheme and back-tracking line search with l-1 norm merit function
% and Armijo condition.
%
%   INPUTS:
%           fun         =   function providing the scalar cost function value 
%                           f(x) and the vectors of p nonlinear equality
%                           constraints g(x), g:R^n->R^p, and q nonlinear 
%                           inequality constraints h(x)
%           x0          =   initial guess for the optimization variables
%           A           =   matrix defining linear equality constraints
%           b           =   vector defining linear equality constraints
%           C           =   matrix defining linear inequality constraints
%           d           =   vector defining linear inequality constraints
%           p           =   number of nonlinear equality constraints
%           q           =   number of nonlinear inequality constraints
%           myoptions   =   optimization options prepared with myoptimset
%                           command
%
%   OUTPUTS:
%           xstar       =   exit value, either a local minimizer or the
%                           value of x at the last iterate
%           fxstar      =   cost function evaluated at xstar
%           niter       =   number of employed iterations 
%           exitflag    =   termination condition:
%                           -2: max iterations reached, unfeasible point
%                           -1: max. number of iterations reached
%                            1: local minimum possible, gradient condition
%                            2: local minimum possible, step size condition
%                            3: local minimum possible, cost decrease condition
%           xsequence   =   sequence of iterations {xk} (only if option
%                           xsequence is set to 'on')

%% Initialization
n               =   length(x0);
k               =   0;
deltaxk_rel     =   1;
deltaf_rel     	=   1;
lambdak         =   zeros(p+size(A,1),1);
muk             =   zeros(q+size(C,1),1);
sigmak          =   zeros(p+size(A,1),1);
tauk            =   zeros(q+size(C,1),1);
eq_constr_max   =   0;
ineq_constr_min =   0;


if ~isempty(myoptions.outputfcn)
    outputfun   =   myoptions.outputfcn;
end
if strcmp(myoptions.display,'Iter')
    fprintf('Iteration       NormGrad          Cost      Equality     Inequality     Rel. cost         Rel. x     Line-search\r')
end

% Start sequence of optimization variables at each iteration
xsequence           =   [];
if strcmp(myoptions.xsequence,'on')
    xsequence       =   [xsequence, x0];
end

%% Iterations
if strcmp(myoptions.Hessmethod,'GN')    % Constrained Gauss-Newton method
    % Compute cost and constraints and their gradients
    funF                =   myoptions.GN_funF;
    xk                  =   x0;
    [Vk,gradVk]         =   mygradient(funF,x0,myoptions.gradmethod,myoptions.graddx);
    Fxk                 =   Vk(1:end-(p+q),1);
    gradFxk             =   gradVk(:,1:end-(p+q));
    if ~isempty(A)
        gxk                 =   [A*xk-b;Vk(end-p-q+1:end-q,1)];
        gradgk              =   [A',gradVk(:,end-p-q+1:end-q,1)];
    else
        gxk                 =   Vk(end-p-q+1:end-q,1);
        gradgk              =   gradVk(:,end-p-q+1:end-q,1);
    end
    if ~isempty(C)
        hxk                 =   [C*xk-d;Vk(end-q+1:end,1)];
        gradhk              =   [C',gradVk(:,end-q+1:end,1)];
    else
        hxk                 =   Vk(end-q+1:end,1);
        gradhk              =   gradVk(:,end-q+1:end,1);
    end
    fxk                 =   Fxk'*Fxk;
    gradfxk             =   2*gradFxk*Fxk;
    
    % Compute graident of the Lagrange function
    gradLagr            =   gradfxk-gradgk*lambdak-gradhk*muk;
    
    % Compute worst-case equality and inequality constraints (for feedback
    % and termination conditions)
    if ~isempty(gxk)
        eq_constr_max       =   max(abs(gxk));
    end
    if ~isempty(hxk)
        ineq_constr_min     =   min(hxk);
    end
    
    % Feedback and output function
    if strcmp(myoptions.display,'Iter')
        if ineq_constr_min<=0 && sign(ineq_constr_min)==-1
            fprintf('%9.0f    %7.5e   %6.5e   %6.5e   %6.5e   %6.5e    %6.5e            %4.0f\r',...
                k,norm(gradLagr),fxk,eq_constr_max,ineq_constr_min,deltaf_rel,deltaxk_rel,0)
        else
            fprintf('%9.0f    %7.5e   %6.5e   %6.5e    %6.5e   %6.5e    %6.5e            %4.0f\r',...
                k,norm(gradLagr),fxk,eq_constr_max,ineq_constr_min,deltaf_rel,deltaxk_rel,0)
        end
    end
    if ~isempty(myoptions.outputfcn)
        outputfun(xk);
    end
    
    % Solver iterations
    while norm(gradLagr) > myoptions.tolgrad...
            && k < myoptions.nitermax...
            && deltaxk_rel > myoptions.tolx...
            && deltaf_rel > myoptions.tolfun...
            || (max(eq_constr_max,-ineq_constr_min) > myoptions.tolconstr)
        
        % Compute search direction
        Hk                  =   2*(gradFxk*gradFxk')+myoptions.GN_sigma*eye(n);
        [pk,~,~,~,LagMult]  =   quadprog(Hk,gradfxk,-gradhk',hxk,gradgk',-gxk,[],[],[],myoptions.QPoptions);
        lambda_tilde        =   -LagMult.eqlin;
        mu_tilde            =   LagMult.ineqlin;
        delta_lambda        =   lambda_tilde-lambdak;
        delta_mu            =   mu_tilde-muk;
        
        % Update merit function
        sigmak              =   max(abs(lambda_tilde),(sigmak+abs(lambda_tilde))/2);
        tauk                =   max(abs(mu_tilde),(tauk+abs(mu_tilde))/2);
        T1fun               =   @(x)MeritT1(funF,x,A,b,C,d,p,q,sigmak,tauk,'GN');
        
        % Compute current merit function value
        T1k                 =   T1fun(xk);

        % compute directional derivative of merit function
        ind_h_viol          =   hxk<=0;
        gradhk_viol         =   gradhk(:,ind_h_viol);
        tauk_viol           =   tauk(ind_h_viol,1);
        DT1k                =   gradfxk'*pk-sigmak'*abs(gxk)-tauk_viol'*gradhk_viol'*pk;
    
        % Line search with merit function and directional derivative
        [xkp1,fxkp1,niter_LS,tk]    =  	linesearch_merit(T1fun,T1k,DT1k,xk,pk,...
            myoptions.ls_tkmax,myoptions.ls_beta,...
            myoptions.ls_c,myoptions.ls_nitermax);
        
        % Update relative changes of cost and optimziation variables
        deltaxk_rel        	=   norm(xkp1-xk)/max(eps,norm(xk));
        deltaf_rel       	=   abs(fxkp1-fxk)/max(eps,abs(fxk));
        
        % Update primal and dual variables
        k                   =   k+1;
        xk                 	=   xkp1;
        lambdak            	=   lambdak+tk*delta_lambda;
        muk                	=   muk+tk*delta_mu;
        
        % Compute new cost, constraints, and their gradients
        [Vk,gradVk]         =   mygradient(funF,xk,myoptions.gradmethod,myoptions.graddx);
        Fxk                 =   Vk(1:end-(p+q),1);
        gradFxk             =   gradVk(:,1:end-(p+q));
        fxk                 =   Fxk'*Fxk;
        gradfxk             =   2*gradFxk*Fxk;
        
        if ~isempty(A)
            gxk                 =   [A*xk-b;Vk(end-p-q+1:end-q,1)];
            gradgk              =   [A',gradVk(:,end-p-q+1:end-q,1)];
        else
            gxk                 =   Vk(end-p-q+1:end-q,1);
            gradgk              =   gradVk(:,end-p-q+1:end-q,1);
        end
        if ~isempty(C)
            hxk                 =   [C*xk-d;Vk(end-q+1:end,1)];
            gradhk              =   [C',gradVk(:,end-q+1:end,1)];
        else
            hxk                 =   Vk(end-q+1:end,1);
            gradhk              =   gradVk(:,end-q+1:end,1);
        end
        
        % Compute gradient of the Lagrangian
        gradLagr            =   gradfxk-gradgk*lambdak-gradhk*muk;
        
        % Compute worst-case equality and inequality constraints (for feedback
        % and termination conditions)
        if ~isempty(gxk)
            eq_constr_max       =   max(abs(gxk));
        end
        if ~isempty(hxk)
            ineq_constr_min     =   min(hxk);
        end
        
        % Feedback and output function
        if strcmp(myoptions.display,'Iter')
        if ineq_constr_min<=0 && sign(ineq_constr_min)==-1
            fprintf('%9.0f    %7.5e   %6.5e   %6.5e   %6.5e   %6.5e    %6.5e            %4.0f\r',...
                k,norm(gradLagr),fxk,eq_constr_max,ineq_constr_min,deltaf_rel,deltaxk_rel,niter_LS)
        else
            fprintf('%9.0f    %7.5e   %6.5e   %6.5e    %6.5e   %6.5e    %6.5e            %4.0f\r',...
                k,norm(gradLagr),fxk,eq_constr_max,ineq_constr_min,deltaf_rel,deltaxk_rel,niter_LS)
        end
        end
        if ~isempty(myoptions.outputfcn)
            outputfun(xk);
        end
        
        % Store sequence of optimization variables at each iteration
        if strcmp(myoptions.xsequence,'on')
            xsequence       =   [xsequence, xk];
        end
    end
    
end

%% Termination
xstar   =   xk;
fxstar  =   fxk;
if max(eq_constr_max,-ineq_constr_min) <= myoptions.tolconstr
    if norm(gradLagr) <= myoptions.tolgrad
        exitflag    =   1;
        if strcmp(myoptions.display,'Iter')
            fprintf('Local minimum possible, directional derivative smaller than tolerance. Constraints satisfied.\r')
        end
    elseif k >= myoptions.nitermax 
        exitflag    =   -1;
        if strcmp(myoptions.display,'Iter')
            fprintf('Maximum number of iterations reached. Constraints satisfied.\r')
        end
    elseif deltaxk_rel <= myoptions.tolx 
        exitflag    =   2;
        if strcmp(myoptions.display,'Iter')
            fprintf('Local minimum possible, relative step size smaller than tolerance. Constraints satisfied.\r')
        end
    elseif deltaf_rel <= myoptions.tolfun 
        exitflag    =   3;
        if strcmp(myoptions.display,'Iter')
            fprintf('Local minimum possible, relative cost decrease smaller than tolerance. Constraints satisfied.\r')
        end
    end
elseif max(eq_constr_max,-ineq_constr_min) > myotions.tolconstr
    if k >= myoptions.nitermax
        exitflag    =   -2;
        if strcmp(myoptions.display,'Iter')
            fprintf('Maximum number of iterations reached. Constraints not satisfied.\r')
        end
    end
end
    