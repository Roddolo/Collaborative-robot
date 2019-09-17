function [q1_dd, q2_dd] = robot_dynamics(tau, state, d)

F_hx = d(1);
F_hz = d(2);
m_load = d(3);
q1  = state(1);
q2 = state(2);
q1_d = state(3);
q2_d = state(4); 

%robot parameters
g   = 9.81;
J1  = 10;   %inertia link 1
J2  = 10;   %inertia link 2
m1  = 50;   %mass link 1
m2  = 50;   %mass link 2
a1  = 1.2;    %lenght link 1
a2  = 1;  %length link 2
ac1  = a1/2; %centre of mass link 1
ac2  = a2/2; %centre of mass link 2

%kinematic functions
% q1          = q(1);
% q2          = q(2);
% q1_d        = q_d(1);
% q2_d        = q_d(2);
c1          = cos(q1);
s1          = sin(q1);
c2          = cos(q2);
s2          = sin(q2);
c12         = cos(q1+q2);
s12         = sin(q1+q2);

J        = [-a1*sin(q1)-a2*sin(q1+q2)  -a2*sin(q1+q2);
             a1*cos(q1)+a2*cos(q1+q2)   a2*cos(q1+q2)];
         
M_tot =  [m1*ac1^2+m2*(a1^2+ac2^2+2*a1*ac2*c2)+J1+J2     m2*(ac2^2+a1*ac2*c2)+J2;
              m2*(ac2^2+a1*ac2*c2)+J2                    m2*ac2^2+J2];
              
% C_tot = [-a1*q2_d*s2*ac2*m2   -a1*q1_d*s2*ac2*m2;
%           -a1*q1_d*s2*ac2*m2             0];- C_tot*[q1_d;q2_d]

         
g_tot = [(m1*ac1+m2*a1)*g*c1+m2*g*ac2*c12;
                 m2*g*ac2*c12];


 %equation of motion of the 2dof manipulator 
q_dd = inv(M_tot)*(tau + J'*[F_hx;F_hz - m_load*g]  - g_tot);
q1_dd = q_dd(1);
q2_dd = q_dd(2);
 
end