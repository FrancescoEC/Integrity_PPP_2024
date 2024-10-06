function plot_latDegree_LonDegree_Hm(timeS, lat, lon, h)

figure
grid on;
set(gcf,'color','white');
subplot 311
plot(timeS,lat,'.-');legend('lat [deg]');
subplot 312
plot(timeS,lon,'.-');legend('lon [deg]');
subplot 313
plot(timeS,h,'.-');legend('hei [m]');
xlabel('sec');