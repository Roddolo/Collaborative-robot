function    [err_vec,q_sim]    =   robot_sim_err(x, q_reference, qd_reference,pos_ref, m_load, Q)

%% Initialize simulation parameters

tau_opt                         =   [x(1,1); x(2,1)];
scaling                         =   size(q_reference,2);         % number of samples 
tau                             =   zeros(5,3);
t_sim                           =   0:0.02:0.08;


for j = 1:5                                                      %1 --> dimension of the sliding window
    tau(j,2:3) = tau_opt';
end
tau(:,1)                        =   t_sim';
tau_sim                         =   tau;
assignin('base', 'tau_sim', tau_sim);                            % assign to workspace matrix disturbances and pass it to simulink

% Simulate the model

sim('Model_sim');
q_sim = [q1_sim q2_sim]';
qd_sim = [q1d_sim q2d_sim]';
F_human = [Fh_x Fh_z]';

%% Simulation of the trajectories

xi_test                         =   zeros(4,5);
for j = 1:5
    [x, z, xd, zd] = direct_kin_test(q_sim(1,j),q_sim(2,j),...
                                    qd_sim(1,j), qd_sim(2,j));         % first 2 rows positions, second 2 speeds
    xi_test(:,j)   = [x, z, xd, zd]';
end
%% Computation of the error to pass to simulink
error_traj                      =   zeros(5,3);
xi_traj                         =   zeros(5,3);
xi_traj(:,2:3)                  =   xi_test(1:2,:)';
xi_traj(:,1)                    =   t_sim';
error_traj(:,1)                 =   t_sim';
error_traj(:,2:3)               =   pos_ref(:,2:3)-xi_traj(:,2:3);
%% Stack errors in one vector

q_ref                           =   [q_reference; qd_reference];
q_test                          =   [q_sim; qd_sim];
err                             =   Q*(q_ref-q_test);
err_vec                         =   [err(1,:) err(2,:) err(3,:) err(4,:)]'/sqrt(scaling);  

%% Assign the parameters to the base workspace

assignin('base', 'xi0', q_sim);                                         
assignin('base', 'xd0', qd_sim);
assignin('base', 'error', err);
assignin('base', 'xi_sim', xi_test);
assignin('base', 'error_xi', error_traj);
assignin('base', 'F_hu', F_human);
end




