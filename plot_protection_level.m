function plot_protection_level(PL_fig,rtk)
global gls motionV2

pos_error = rtk.sol.pos - [motionV2.Pos_X(gls.ti+1),motionV2.Pos_Y(gls.ti+1),motionV2.Pos_Z(gls.ti+1)];


figure(PL_fig)
for i = 1:3
subplot(3,1,i)
hold on
plot(gls.ti,gls.Log.PPP.PL(i),'b*');
plot(gls.ti,pos_error(i),'ro')
end

end