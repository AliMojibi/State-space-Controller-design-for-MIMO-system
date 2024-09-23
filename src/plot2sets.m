function plot2sets(t, x, x_lin)
figure
subplot(321);plot(t, x(:, 1),'-b ', 'LineWidth',2);hold on; plot(t, x_lin(:, 1), '--r ', 'LineWidth',2);xlabel('t(s)');ylabel('x(m)');grid on;legend('Non-linear sys.', 'Linear sys.')
subplot(322);plot(t, x(:, 2),'-b ', 'LineWidth',2);hold on; plot(t, x_lin(:, 2), '--r ', 'LineWidth',2);xlabel('t(s)');ylabel('x dot(m)');grid on;legend('Non-linear sys.', 'Linear sys.')
subplot(323);plot(t, x(:, 3),'-b ', 'LineWidth',2);hold on; plot(t, x_lin(:, 3), '--r ', 'LineWidth',2);xlabel('t(s)');ylabel('theta1(rad)');grid on;legend('Non-linear sys.', 'Linear sys.')
subplot(324);plot(t, x(:, 4),'-b ', 'LineWidth',2);hold on; plot(t, x_lin(:, 4), '--r ', 'LineWidth',2);xlabel('t(s)');ylabel('theta1 dot(rad/s)');grid on;legend('Non-linear sys.', 'Linear sys.')
subplot(325);plot(t, x(:, 5),'-b ', 'LineWidth',2);hold on; plot(t, x_lin(:, 5), '--r ', 'LineWidth',2);xlabel('t(s)');ylabel('theta2(rad)');grid on;legend('Non-linear sys.', 'Linear sys.')
subplot(326);plot(t, x(:, 6),'-b ', 'LineWidth',2);hold on; plot(t, x_lin(:, 6), '--r ', 'LineWidth',2);xlabel('t(s)');ylabel('theta2 dot(rad/s)');grid on;legend('Non-linear sys.', 'Linear sys.')
set(findobj(gcf,'type','axes'),'FontName','TimesNewRoman','FontSize',14,'FontWeight','Bold', 'LineWidth', 1,'layer','top');
end