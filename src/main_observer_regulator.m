%----------------------------------------
% Ali Mojibi 401202618                   |
% Advanced Automatic Control             |
% Double invereted Pendulum Homework     |
% 12/10/2023                             |
%----------------------------------------
clc, clear all, close all

g = 9.81;
m_0 = 2; m = 0.5; L = 0.5; % SI units

% ----------------------------- part 03 -----------------------------------
A = A_lin(L,g,m,m_0)
B = B_lin(L,m,m_0)
pause
mm = rank(B) % num of independent columns in B
nn = length(A) % num of state variables
% n - m = 4
M_c = [B, A*B, A^2*B, A^3*B, A^4*B]% Controlability matrix
rank(M_c) % M_c is full rank --> OK

% ----------------------------- part 04 -----------------------------------
% y = [x theta_1]
C_1 = [1 0 0 0 0 0;0 0 1 0 0 0];
N_1 = [C_1; C_1*A; C_1*A^2; C_1*A^3; C_1*A^4]
rank(N_1) % full rank --> Observable
C = C_1;
% ----------------------------- part 05 -----------------------------------
% y = [theta_1 theta_2]
C_2 = [0 0 1 0 0 0;0 0 0 0 1 0];
N_2 = [C_2; C_2*A; C_2*A^2; C_2*A^3; C_2*A^4]
rank(N_2) %rank=4 --> Not Observable
% ----------------------------- part 06 -----------------------------------
% EESA
mu(1) = -2; mu(2) = -2; mu(3) = -2-j; mu(4) = -2+j; mu(5) = -3; mu(6) = -3;
S1 = null([A-mu(1)*eye(6) B])
S3 = null([A-mu(3)*eye(6) B])
S4 = null([A-mu(4)*eye(6) B])
S5 = null([A-mu(5)*eye(6) B])
v1 = S1(1:6, 1); q1 = S1(7:end, 1);
v2 = S1(1:6, 2); q2 = S1(7:end, 2);
v3 = S3(1:6, 1); q3 = S3(7:end, 1);
v4 = S4(1:6, 1); q4 = S4(7:end, 1);
v5 = S5(1:6, 1); q5 = S5(7:end, 1);
v6 = S5(1:6, 2); q6 = S5(7:end, 2);

K_EESA = -[q1 q2 q3 q4 q5 q6]*inv([v1 v2 v3 v4 v5 v6])
K_EESA=real(K_EESA);
eig(A-B*K_EESA)
% ----------------------------- part 07 -----------------------------------
% Observer Design using GCCf algorithm
A = A.'; B=C_1.';
mu = -10*ones(1, 6);
M_h = [B(:, 1) A*B(:, 1) A^2*B(:, 1) A^3*B(:, 1), B(:, 2) A*B(:, 2)] % gamma1 = 4, gamma2 = 2
M_hi = inv(M_h);
e14T=M_hi(4, :)
e22T=M_hi(4+2, :)
Ti = [e14T; e14T*A; e14T*A^2; e14T*A^3; e22T; e22T*A]
T=inv(Ti)
% dX = AX + Bu --> dZ = A_G Z + B_G v; X=TZ u=fw w = H - vZ 
Ag1 = [zeros(3, 1), eye(3);zeros(1, 4)];
Ag2 = [0 1;0 0];
A_G = [Ag1, zeros(4, 2);zeros(2, 4) Ag2]
Bg1 = [0 0 0 1].';
Bg2 = [0 1].';
B_G = [Bg1 zeros(4,1);zeros(2,1) Bg2];
% Ti*B*F = B_G --> B*F = T*B_G-->(B.'*B)*F=B.'*T*B_G-->F = inv(B.'*B)*B.'*T*B_G
F = inv(B.'*B)*B.'*T*B_G;
H = B_G.'*(A_G-Ti*A*T);

%desired polynomials
syms s
P1 = flip(coeffs(expand((s-mu(1))^4)));% gamma1 = 4
P2 = flip(coeffs(expand((s-mu(1))^2)));% gamma2 = 2
Ad1 = [zeros(3,1), eye(3);-flip(P1(2:end))]
Ad2 = [0 1;-flip(P2(2:end))]
A_d =[Ad1, zeros(4,2);zeros(2, 4), Ad2]

Gamma = B_G.'*(A_G-A_d)
K_o = double(F*(Gamma-H)*Ti).'

filename = 'K_o.mat';
% Save the matrix to the file
save(filename, 'K_o');

% ----------------------------- part 08 -----------------------------------
A = A_lin(L,g,m,m_0);
B = B_lin(L,m,m_0);
eig(A-B*K_EESA)
t_span = [0:.01:10];
% non linear system
ode_fun = @(t, x)system_dynamics(t, x, K_EESA,L,g,m,m_0);
x0 = [0 0.1 5*pi/180 3*pi/180 -2*pi/180 -5*pi/180]
[t, x] = ode45(ode_fun, t_span, x0);

%linear system
ode_fun = @(t, x)system_dynamics_lin(t, x, K_EESA,L,g,m,m_0);
x0 = [0 0.1 5*pi/180 3*pi/180 -2*pi/180 -5*pi/180]
[t, x_lin] = ode45(ode_fun, t_span, x0);


plot2sets(t, x, x_lin)
% ----------------------------- part 09 -----------------------------------
t_span = [0:.01:60];
%linear system
ode_funl = @(t, x)observed_system_dynamics_lin(t, x, K_EESA,K_o,C,L,g,m,m_0, 0);
x0 = [0 0.1 5*pi/180 3*pi/180 -2*pi/180 -5*pi/180, 0 0.1 5*pi/180 3*pi/180 -2*pi/180 -5*pi/180];
[t, x_lin] = ode45(ode_funl, t_span, x0);

% non linear system
ode_fun = @(t, x)observed_system_dynamics(t, x, K_EESA,K_o,C,L,g,m,m_0, 0);
x0 = [0 0.1 5*pi/180 3*pi/180 -2*pi/180 -5*pi/180, 0 0.1 5*pi/180 3*pi/180 -2*pi/180 -5*pi/180];
[t, x] = ode45(ode_fun, t_span, x0);
plot2sets(t, x, x_lin)



