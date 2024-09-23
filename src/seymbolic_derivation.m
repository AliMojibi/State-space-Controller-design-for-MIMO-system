clc 
clear all

syms t real
syms x theta_1 theta_2 real
syms x_dot theta_1_dot theta_2_dot real
%syms q1 q2 q3 real % Generalized Coordinates q1 = x, q2 = theta_1, q3 = theta_2 
%syms dq1 dq2 dq3 real% Generalized Velocity dq/dt
syms m_0 m L g real real % System Parameters
syms u tau real% External Forces and Torque
q1 = x; q2 = theta_1; q3 = theta_2;
dq1 = x_dot; dq2 = theta_1_dot; dq3 = theta_2_dot;

I_1 = 1/12*m*L^2;
I_2 = I_1;
%% Velocity Calculation
v_A = [dq1; 0; 0];
omega_1 = dq2*[0; 0; 1];
r_G1A = L/2*[sin(q2), -cos(q2), 0].';
v_G1 = v_A + cross(omega_1, r_G1A)
r_BA = L*[sin(q2), -cos(q2), 0].';
v_B = v_A + cross(omega_1, r_BA);
omega_2 = dq3*[0; 0; 1];
r_G2B = L/2*[sin(q3), -cos(q3), 0].';
v_G2 = v_B + cross(omega_2, r_G2B)

%% Kinetic Energy
T_cart = 1/2*m_0*v_A.'*v_A;
T_1 = 1/2*m*v_G1.'*v_G1 + 1/2*I_1*omega_1.'*omega_1;
T_2 = 1/2*m*v_G2.'*v_G2 + 1/2*I_2*omega_2.'*omega_2;

T = T_cart + T_1 + T_2


%% Potential Energy
r_A = q1*[1 0 0].';
r_G1 = r_A + r_G1A;
r_G2 = r_A + r_BA + r_G2B

V_G = -[0 1 0]*(m*g*r_G1+m*g*r_G2); % we need the y component of the r_Gi

V = V_G;

% Lagrangian
L = T - V;


q = [q1 q2 q3].';
dq = [dq1 dq2 dq3].';

% Determination of Generalized Forces
F1 = u*[1 0 0].';
Tau1 = tau*[0 0 1].';


Power = F1.'*v_A + Tau1.'*(omega_1-omega_2);

Q = jacobian(Power,dq)';

dL_dq = jacobian(L,q);  % dL/dq
dL_ddq = jacobian(L,dq); % dL/d(qdot)

M = jacobian (dL_ddq,dq);
Bias = jacobian(dL_ddq,q)*dq + diff(dL_ddq,t)' - dL_dq';

M = simplify(M);
Bias = simplify(Bias);
Q = simplify(Q);

ddq = simplify(M\(Q-Bias))

% --------------------------State space Modeling---------------------------
X = [x x_dot theta_1 theta_1_dot theta_2 theta_2_dot]% state variables
X_dot = [x_dot, ddq(1), theta_1_dot, ddq(2), theta_2_dot, ddq(3)].'; % deriviative of state variables
u_vec = [u tau].';% inputs vectors
B = simplify(jacobian(X_dot, u_vec))
A = simplify(jacobian(X_dot, X))

A_lin = subs(A, [x x_dot theta_1 theta_1_dot theta_2 theta_2_dot], zeros(1, 6))
B_lin = subs(B, [x x_dot theta_1 theta_1_dot theta_2 theta_2_dot], zeros(1, 6))

% matlabFunction(A, "file", "A_non_lin")
% matlabFunction(B, "file", "B_non_lin")
% matlabFunction(A_lin, "file", "A_lin")
% matlabFunction(B_lin, "file", "B_lin")
% 
