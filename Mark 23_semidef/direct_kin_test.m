function [xe, ze, x_d, z_d]   = direct_kin_test(q1, q2, q1_d, q2_d)
%robot parameters
a1  = 1.2;    %lenght link 1
a2  = 1;  %length link 2

%kinematic functions
c1= cos(q1);
s1= sin(q1);
c2= cos(q2);
s2= sin(q2);
c12= cos(q1+q2);
s12= sin(q1+q2);

% Computation of end effector quantities given the measure joints
% quantities
if q1 == 0 && q2 == 0
    xe = 0;
    ze = 0;
   
else 
    xe = 0.5 + a1*c1 + a2*c12;
    ze = a1*s1 + a2*s12;
    
end
%xi = [xe;ze];
    
%% computations of speeds
J        = [-a1*sin(q1)-a2*sin(q1+q2)  -a2*sin(q1+q2);
             a1*cos(q1)+a2*cos(q1+q2)   a2*cos(q1+q2)];

xi_d = J*[q1_d;q2_d];

x_d = xi_d(1,1);
z_d = xi_d(2,1);
   
end
