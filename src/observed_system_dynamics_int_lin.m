function dx = observed_system_dynamics_int_lin(t, x, K,K_o,C,L,g,m,m_0)
    dx = zeros(size(x));
    r=x(1); r_dot=x(2);theta_1 = x(3); theta_1_dot=x(4);theta_2 = x(5); theta_2_dot = x(6);
    xi_1 = x(7); xi_2 = x(8);
    r_h=x(9); r_dot_h=x(10);theta_1_h = x(11); theta_1_dot_h=x(12);theta_2_h = x(13); theta_2_dot_h = x(14);
    
    X = x(1:6);X_h = x(9:end);
    
    
    ref_r = sign(sin(0.2*t));
    ref_theta_1 = 5*pi*sign(sin(0.2*t))/180;
    ref = [ref_r;ref_theta_1];

    U = -K*[X_h; xi_1; xi_2];
    u = U(1); tau=U(2);

    A = A_lin(L,g,m,m_0);
    B = B_lin(L,m,m_0);

    A_h = A;
    B_h = B;
    
    dX = A*X + B*U;
    xi_dot = ref-C*X;
    dX_h = A_h * X_h + B_h*U + K_o*C*(X-X_h);

    dx = [dX; xi_dot; dX_h]
    

end