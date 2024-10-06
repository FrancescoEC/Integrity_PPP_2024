function plot_ENU(timeS, E, N, U )

figure
set(gcf,'Color','white');
grid on;
subplot 311
plot(timeS, E,'.-'); legend('E [m]')
subplot 312
plot(timeS, N,'.-'); legend('N [m]')
subplot 313
plot(timeS, U,'.-'); legend('U [m]')

figure
set(gcf,'Color','white');
grid on;
plot(E,N,'.-k');
xlabel('E [m]');
ylabel('N [m]');
grid on;