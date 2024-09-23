function dx = system_dynamics_lin(t, x, K,L,g,m,m_0)
    dx = zeros(size(x));
    %theta_1 = x(3); theta_1_dot=x(4);theta_2 = x(5); theta_2_dot = x(6)
    U = -K*x;
    u = U(1); tau=U(2);

    A = A_lin(L,g,m,m_0);
    B = B_lin(L,m,m_0);

    dx = A*x+B*[u;tau];
    

end