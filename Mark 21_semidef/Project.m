% Constrained numerical optimization for estimation and control. 
% Project by Rodrigo Senofieni, Andrea Ghezzi and Luca Brolese. 
% AA 2018/2019
% Real time estimation of the human force acting on the end-effector of a collaborative robot. 
% Estimation procedure is carried out by the usage of a MHE approach based on past records of 
% the error (between reference and simulated joints position), evaluated in a sliding window of dimension M.

%%
close all 
clear all 
clc

%% Parameters
fs                              =   50;                         % sampling frequency 
Ts                              =   1/fs;                       % sampling period
tau_p                           =   0.1;                        % pole of the human TF
tau_z                           =   1;                          % zero of the human TF
k                               =   100;                        % Gain of the human TF
m_load                          =   10;                         % Mass of the load
Fmax                            =   200;                        % Maximum force of the human. Can be used in saturation block
s                               =   tf('s');                    % Initialize the TF
H                               =   k*(1+tau_z*s)/(1+tau_p*s);  % Human TF
[numH, denH]                    =   tfdata(H, 'v');             % Human TF parameters for simulink
% Parameters of the robot
g                               =   9.81;                       % gravitational acceleration
J1                              =   10;                         % inertia link 1
J2                              =   10;                         % inertia link 2
m1                              =   50;                         % mass link 1
m2                              =   50;                         % mass link 2
a1                              =   1.2;                        % lenght link 1
a2                              =   1;                          % lenght link 2
l1                              =   a1/2;                       % centre of mass link 1
l2                              =   a2/2;                       % centre of mass link 2

%% Data acquisition

sim('Model_daq');
M                               =   5;                                         % dimension of the sliding window
downsampling                    =   50;                                        % Reduce sampling frequency for higher simulation speed with FFD                 
q_reference                     =   [state(1:Ts*downsampling:500,1)';          % q1
                                    state(1:Ts*downsampling:500,2)'];          % q2
qd_reference                    =   [state(1:Ts*downsampling:500,3)';          % q1d
                                    state(1:Ts*downsampling:500,4)'];          % q2d                
Fh_reference                    =   [F_h(1:Ts*downsampling:500,1)';            % Fh_x_ref
                                     F_h(1:Ts*downsampling:500,2)'];           % Fh_z_ref
position_ref                    =   [xi(1:Ts*downsampling:500,1)';             % x_ref
                                     xi(1:Ts*downsampling:500,2)'];            % z_ref
speed_ref                       =   [xid(1:Ts*downsampling:500,1)';            % xd_ref
                                     xid(1:Ts*downsampling:500,2)'];           % zd_ref                   

%% Parameters for the optimization routine
Q                               =    eye(4);               % weight matrix diag([1 0.8 0.6 0.4 0.2]);
x0                              =    [0;-200];             % initialize parameter estimate [F_hx; F_hz]
scaling                         =    Ts*downsampling;      % set up the scaling for the downsampling
dim                             =    500/scaling;
F_h_star                        =    zeros(2,dim);         % initialize vector of F_h_estimated
tau                             =    zeros(2,M);           % dimension of the torque vector passed to the functions. Torques are setted to zero
xi_test                         =    zeros(4,500);         % vector for storing the simulated trajectories
xi0                             =    zeros(2,M);
xi0(:,2)                        =    [2.18;2.72];          % Initial condition for the integrators, position
xd0                             =    zeros(2,M);           % Initial condition for the integrators, speed

%% Optimization Routine
% Linear equality constraint parameters
A                               =   [];
b                               =   [];
% Linear inequality constraint parameters
C                               =   [];                    %1 0; -1 0;0 1; 0 -1
d                               =   [];                    %-5000;-5000;-5000;-5000

for i = M+1:dim
%Initialize solver options
myoptions                       =   myoptimset;
myoptions.Hessmethod            =	'GN';
myoptions.GN_funF               =	@(x)robot_sim_err(x,q_reference(:,i-M:i-1),...
                                                      qd_reference(:,i-M:i-1), tau, m_load, Q);
myoptions.gradmethod            =	'CD';
myoptions.graddx                =	2^-15;
myoptions.tolgrad               =	1e-7;
myoptions.ls_nitermax           =	50;
myoptions.xsequence             =	'on';

%initial condition for integrators
q1_in                           =   xi0(1,2)
q2_in                           =   xi0(2,2)
q1_d_in                         =   xd0(1,2);
q2_d_in                         =   xd0(2,2);

% Run solver
[xstar,fxstar,niter,exitflag]   =   myfmincon(@(x)robot_sim_cost(x,q_reference(:,i-M:i-1),...
                                            qd_reference(:,i-M:i-1),tau, m_load, Q),...
                                            x0,A,b,C,d,0,0,myoptions);
xstar

% Store the obtained result of xstar
F_h_star(:,i)                   =   xstar;
xi_test(:,i)                    =   xi_sim(:,1);

%xi0 variable assigned from the function robot_sim_err, 
iter = i
error

end

%% Plot the results

t          =        0:0.02:9.98;
figure
subplot 221 
plot (t, xi_test(1,:)),grid on
hold on
plot (t,position_ref(1,:)),grid on
legend('simulated x trajectory','reference x trajectory')

subplot 222 
plot (t,xi_test(2,:)),grid on
hold on
plot (t,position_ref(2,:)),grid on
legend('simulated z trajectory','reference z trajectory')


subplot 223
plot (t,xi_test(3,:)),grid on
hold on
plot (t,speed_ref(1,:)), grid on
legend('simulated x speed','reference x speed')

subplot 224
plot (t,xi_test(4,:)),grid on
hold on
plot (t,speed_ref(2,:)), grid on
legend('simulated z speed','reference z speed')

figure
subplot 211
plot (t, F_h_star(1,:)), grid on
hold on
plot (t, Fh_reference(1,:)), grid on
legend('simulated F_hx','reference Fh_x')

subplot 212
plot (t, F_h_star(2,:)), grid on
hold on
plot (t, Fh_reference(2,:)), grid on
legend('simulated Fh_z','reference Fh_z')














