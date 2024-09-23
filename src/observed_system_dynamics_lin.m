function dx = observed_system_dynamics_lin(t, x, K,K_o,C,L,g,m,m_0, u_f)
    dx = zeros(size(x));
    r=x(1); r_dot=x(2);theta_1 = x(3); theta_1_dot=x(4);theta_2 = x(5); theta_2_dot = x(6);
    r_h=x(7); r_dot_h=x(8);theta_1_h = x(9); theta_1_dot_h=x(10);theta_2_h = x(11); theta_2_dot_h = x(12);
    
    X = x(1:6); X_h = x(7:end);
    
    U_ff = zeros(2,1);
    
    if u_f
      ref_r = sign(sin(0.2*t));
      ref_theta_1 = 5*pi*sign(sin(0.2*t))/180;

      A = A_lin(L,g,m,m_0);
      B = B_lin(L,m,m_0);
      U_ff = -inv(C*inv(A-B*K)*B)*[ref_r;ref_theta_1]; % feed forward value
    end

    U = -K*X_h+U_ff;


    A_l = A_lin(L,g,m,m_0);
    B_l = B_lin(L,m,m_0);

    A_h = A_l;
    B_h = B_l;

    dX = A_l*X + B_l*U;
    dX_h = A_h*X_h + B_h*U + K_o*C*(X-X_h);

    dx = [dX; dX_h]
    

end