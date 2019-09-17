function    [err_vec,q_sim]    =   robot_sim_err(x, q_reference,qd_reference,tau, m_load, Q)

%% Initialize simulation parameters

F_h                             =   [x(1,1); x(2,1)];            % set the optimization variables as [Fh_x;Fh_z]
scaling                         =   size(q_reference,2);         % number of samples 
disturbances                    =   zeros(5,4);
t_sim                           =   0:0.02:0.08;

for j = 1:5                                                      %1 --> dimension of the sliding window
    disturbances(j,2:3) = F_h';
    disturbances(j,4) = m_load;
end
disturbances(:,1)               =   t_sim';
dist_sim                        =   disturbances;                % assign to workspace matrix disturbances and pass it to simulink
assignin('base', 'dist_sim', dist_sim);

tau0                            =   [t_sim' tau'];
tau_sim                         =   tau0;
assignin('base', 'tau_sim', tau_sim);                            % assign to workspace matrix disturbances and pass it to simulink

% Simulate the model

sim('Model_sim');
q_sim = [q1_sim q2_sim]';
qd_sim = [q1d_sim q2d_sim]';

%% Simulation of the trajectories

xi_test                         =   zeros(4,5);
for j = 1:5
    [x, z, xd, zd] = direct_kin_test(q_sim(1,j),q_sim(2,j),...
                                    qd_sim(1,j), qd_sim(2,j));         % first 2 rows positions, second 2 speeds
    xi_test(:,j)   = [x, z, xd, zd]';
end

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
end




