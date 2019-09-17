function [T1val] = MeritT1(fun,xk,A,b,C,d,p,q,sigma,tau,mode)
% MERITT1 Computes the vlaue of the merit function.
%
%   INPUTS:
%           fun         =   cost function
%           xk          =   current point
%           A           =   matrix defining linear equality constraints
%           b           =   vector defining linear equality constraints
%           C           =   matrix defining linear inequality constraints
%           d           =   vector defining linear inequality constraints
%           p           =   number of nonlinear equality constraints
%           q           =   number of nonlinear inequality constraints
%           sigma       =   vector of weights, one for each equality
%                           constraint
%           tau         =   vector of weights, one for each inequality
%                           constraint
%           mode        =   'GN': Gauss-Newton
%                           'BFGS': BFGS  
%   OUTPUTS:
%           T1val       =   value of the merit function

[Vk]	=   fun(xk);
if strcmp(mode,'GN')
    Fxk     =   Vk(1:end-(p+q),1);
    fxk   	=   Fxk'*Fxk;
elseif strcmp(mode,'BFGS')
    fxk     =   Vk(1:end-(p+q),1);
end
if ~isempty(A)
    gxk                 =   [A*xk-b;Vk(end-p-q+1:end-q,1)];
else
    gxk                 =   Vk(end-p-q+1:end-q,1);
end
if ~isempty(C)
    hxk                 =   [C*xk-d;Vk(end-q+1:end,1)];
else
    hxk                 =   Vk(end-q+1:end,1);
end
T1val  	=   fxk+sigma'*abs(gxk)+tau'*max(zeros(size(hxk)),-hxk);