sp = skyplot([],[]);
for tIdx = 1:LEO_PNT(1).numTimeStamps
    set(sp, ...
        AzimuthData=LEO_PNT_T.azTotal(tIdx,:), ...
        ElevationData=LEO_PNT_T.elTotal(tIdx,:), ...
        LabelData=LEO_PNT_T.NameTotal,ColorOrder=(tIdx*40/256)*ones(1,3));
    % For slower animation, use 'drawnow' instead of 'drawnow limitrate'
    drawnow 
    pause(0.5)
end

% Find the PRN index of each satellite
satNames = char(sat(:).Name');
prnIndex = double(string(satNames(:,5:end)));

% To better visualize each GPS satellite, scale the status with the PRN
% index
acStatsAllTime = double(acStatsAllTime);
acStatsAllTime(acStatsAllTime == 0) = NaN;

% Plot the satellite visibility chart
colors = colororder;
firstLineColor = colors(1,:);
plot(timeHistory,prnIndex.*acStatsAllTime, ...
    Color=firstLineColor,LineWidth=1)
xlim([timeHistory(1) timeHistory(end)])
ylim([min(prnIndex)-1 max(prnIndex)+1])
xlabel("Time")
ylabel("Satellite PRN Index")
title("Satellite Visibility Chart")
yticks(prnIndex)
grid on
