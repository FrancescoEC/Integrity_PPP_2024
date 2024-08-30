function y=periodic_parabola(a0, a1, a2,phase, T, totalT)
    % Check if inputs are valid
    if T <= 0
        error('Period T must be greater than zero.');
    end
    if totalT <= 0
        error('Total time totalT must be greater than zero.');
    end
    if totalT < T
        warning('Total time totalT is less than one period T. Output will be truncated.');
    end

    % Time vector with a reasonable number of points per period
    dt = 1; % You can adjust the resolution by changing the divisor
    t = 0:dt:totalT-1;

    % Function to calculate the parabola within one period
    function y = parabola_within_period(t)
        Index=floor(t/T)+1;
        % This makes sure that we are only considering the time within one period
        t = mod(t, T);
        % Calculate the parabolic function
        y = a0(Index) + a1(Index) * t + 1/2*a2(Index) * t.^2+a1(t+1)+a2(t+1)+a0(t+1)/10*sin(2*pi*(t+1+phase(Index)*T)/T);
        y =abs(y);
    end

    % Apply the parabola function to the entire time vector
    y = arrayfun(@parabola_within_period, t);

%     % Plot the result
%     figure;
%     plot(t, y);
%     title('Periodically Repeated Positive Parabola');
%     xlabel('Time');
%     ylabel('Amplitude');
%     grid on;
end
