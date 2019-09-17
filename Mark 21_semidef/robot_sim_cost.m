function [cost, q_sim] = robot_sim_cost(x,q_meas,qd_ref, m_load,Q)

%Compute the quadratic cost of the error obteined after the simulation
[err_vec, q_sim]    =   robot_sim_err(x, q_meas,qd_ref, m_load, Q);
cost                =   sum(err_vec.*err_vec);

end