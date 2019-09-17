% Constrained numerical optimization for estimation and control. 
% Project by Rodrigo Senofieni, Andrea Ghezzi and Luca Brolese. 
% AA 2018/2019
% Real time optimal control of the torques of a robot, helping the human lifting a load 
% optimal control problem is carried out 

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
Fmax                            =   100;                        % Maximum force of the human. Can be used in saturation block
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
q_reference                     =   [q_ref(1:Ts*downsampling:500,1)';          % q1
                                    q_ref(1:Ts*downsampling:500,2)'];          % q2
qd_reference                    =   [qd_ref(1:Ts*downsampling:500,1)';         % q1d
                                    qd_ref(1:Ts*downsampling:500,2)'];         % q2d                
position_ref                    =   [xi_ref(1:Ts*downsampling:500,1)';             % x_ref
                                     xi_ref(1:Ts*downsampling:500,2)'];            % z_ref
speed_ref                       =   [xid_ref(1:Ts*downsampling:500,1)';            % xd_ref
                                     xid_ref(1:Ts*downsampling:500,2)'];           % zd_ref                   

%% Parameters for the optimization routine
Q                               =    eye(4);               % weight matrix diag([1 0.8 0.6 0.4 0.2]);
x0                              =    [100;100];            % initialize parameter estimate [F_hx; F_hz; tau1; tau2]
scaling                         =    Ts*downsampling;      % set up the scaling for the downsampling
dim                             =    500/scaling;        
F_h_star                        =    zeros(2,dim);         % initialize vector of F_h_estimated
xi_test                         =    zeros(4,500);         % vector for storing the simulated trajectories
xi0                             =    zeros(2,M);
xi0(:,2)                        =    [2.18;2.72];          % Initial condition for the integrators, position
xd0                             =    zeros(2,M);           % Initial condition for the integrators, speed
pos_ref                         =    zeros(5,3);           % initialize empty vector for storing reference position
Fh                              =    zeros(2,500);         % empty vector for storing the resulting human force
sim_time                        =    0:0.02:0.08;          % vector for simulation time
error_xi                        =    zeros(5,3);                  %error along x-z direction for the human closed loop
error_xi(:,2:3)                 =    position_ref(:,1:5)';        %initialize the error = pos_ref
error_xi(:,1)                   =    sim_time';
%% Optimization Routine
% Linear equality constraint parameters
A                               =   [];
b                               =   [];
% Linear inequality constraint parameters
C                               =   [];            
d                               =   [];                  
tic
for i = 1:dim-M
%Initialize solver options
myoptions                       =   myoptimset;
myoptions.Hessmethod            =	'GN';
myoptions.GN_funF               =	@(x)robot_sim_err(x,q_reference(:,i+1:M+i),...
                                                      qd_reference(:,i+1:M+i),pos_ref, m_load, Q);
myoptions.gradmethod            =	'CD';
myoptions.graddx                =	2^-13;
myoptions.tolgrad               =	1e-6;
myoptions.tolx                  =	1e-10;       % Termination tolerance on the relative
                                                 % change of optimization variable
myoptions.ls_nitermax           =	1e2;
myoptions.xsequence             =	'on';

%initial condition for integrators
q1_in                           =   xi0(1,2);
q2_in                           =   xi0(2,2);
q1_d_in                         =   xd0(1,2);
q2_d_in                         =   xd0(2,2);

%initialize at each cycle the position reference for simulink

pos_ref(:,2:3)                  =   position_ref(:,i+1:M+i)';
pos_ref(:,1)                    =   sim_time';

% Run solver
[xstar,fxstar,niter,exitflag]   =   myfmincon(@(x)robot_sim_cost(x,q_reference(:,i+1:M+i),...
                                            qd_reference(:,i+1:M+i),pos_ref, m_load, Q),...
                                            x0,A,b,C,d,0,0,myoptions);
xstar

% Store the obtained result of data
xi_test(:,i)                    =   xi_sim(:,1);
tau_star(:,i)                   =   xstar(1:2,1);
Fh(:,i)                         =   F_hu(:,1);

%xi0 variable assigned from the function robot_sim_err, 
iter = i
error

end
toc
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
plot (t, Fh(1,:)), grid on
hold on
% plot (t, Fh_reference(1,:)), grid on
% legend('simulated F_hx','reference Fh_x')

subplot 212
plot (t, Fh(2,:)), grid on
hold on
% plot (t, Fh_reference(2,:)), grid on
% legend('simulated Fh_z','reference Fh_z')

figure
subplot 211
plot (t(:,1:216), tau_star(1,:)), grid on           %SET THE DIMENSION OF THE time vector, depending on the dimension of tau_star
subplot 212
plot (t(:,1:216), tau_star(2,:)), grid on















