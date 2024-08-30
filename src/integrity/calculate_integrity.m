function [PL1,PL2,PL3,PL4,PRN_list] = calculate_integrity(X,H,P,R,v,Buffer_V_inv,Buffer_y,PRN,PRN_list)
% Created by Heiko Engwerda (NLR)
% 01/07/2024
%
% Function that determines protection levels based on two integrity
% algorithms: AIME (sliding window residual based) and ARAIM (solution
% separation)

% Set constants
PHMI = 1e-7; % probability of misleading hazardous information (per hour)
P_fa = 1e-5; % probability of false alarm
P_THRES = 0.9*PHMI; % threshold on not-monitored faults
minimum_detectable_bias = 0.5; % minimum detectable bias (m)
p_fault = 1e-4; % satellite fault probability (per hour)
p_sat = (1/300)/(30*24); % assuming one satellite failure per month with 300 satellites in service
% p_sat = 1e-5; % satellite fault probability
p_const = 1/(365*24); % assuming one constellation wide failure per year
% p_const = 1e-4; % constellation wide fault probability
epsilon = 1e-6; % Regularization parameter

% Fault detection and exclusion
% obtain innovation covariance
S = H * P * H' + R;

% obtain Kalman gain
K = P * H'* inv(S);

% Sliding window test-stat (AIME)
% compare sliding window PRN with current PRN
[~,idx_rem_sats,idx_add_sats] = setxor(PRN_list,PRN);

% remove satellites from buffer if needed
if ~isempty(idx_rem_sats)
    Buffer_V_inv(idx_rem_sats,idx_rem_sats,:) = [];
    Buffer_y(idx_rem_sats,:) = [];
    PRN_list(idx_rem_sats) = [];
end

% add new satellites to the list
if ~isempty(idx_add_sats)
    PRN_list = [PRN_list,PRN(idx_add_sats)];
end

% determine the order of the satellites in the measurement
[~,idx_meas,idx_buff] = intersect(PRN,PRN_list);

% determine number of satellites
numsat = numel(PRN_list);

Buffer_V_inv(idx_buff,idx_buff,1) = inv(S(idx_meas,idx_meas)+epsilon* eye(size(S)));
Buffer_y(idx_buff,1) = v(idx_meas);

V_avg_inv = zeros(numsat,numsat);
y_avg = zeros(numsat,1);
for window = 1:size(Buffer_V_inv,3)
    V_avg_inv = V_avg_inv + Buffer_V_inv(idx_buff,idx_buff,window);
    y_avg = y_avg + Buffer_V_inv(idx_buff,idx_buff,window)*Buffer_y(idx_buff,window);
end

% detection statistic
s_avg = sqrt(y_avg'*V_avg_inv*y_avg);

% detection threshold
DOF = numsat;
threshold = chi2inv(1-P_fa,DOF);

PL1 = zeros(1,3);
if threshold>0
    if s_avg^2>threshold
        warning('fault detected')
        PL1 = [NaN NaN NaN]';
        % end
    else

        slope_max = sqrt(max(eig(K'*K*S)));

        k = sqrt(2) * erfinv(1-PHMI);

        PL1 = slope_max * sqrt(threshold) + k * sqrt(diag(P(1:3,1:3))); !TODO: Heiko check P indices (should match position)
    end
else
    PL1 = [NaN NaN NaN]';
end

% Single epoch test-stat
test_stat = v' * inv(S) * v;

% Failure detection
if test_stat>threshold
    warning('T>tresh: FDE attempt')
end

% weighting matrix
W = inv(R);

% worst satellite fault vector
[~,worst_sat] = max(diag(W));
fault_vec = zeros(numsat,1);
fault_vec(worst_sat) = minimum_detectable_bias;

% noncentrality parameter
lambda = fault_vec'*W*fault_vec;

% Slope calculation
valid_col = find(abs(sum(H,1))>0);
A = inv(H(:,valid_col)'*W*H(:,valid_col))*H(:,valid_col)'*W;
S = eye(numsat) - H(:,valid_col)*A;

slope_1 = zeros(numsat,1); slope_2 = zeros(numsat,1); slope_3 = zeros(numsat,1);
for sat=1:numsat
    slope = sqrt(A(1,sat)^2+A(2,sat)^2)/sqrt(S(sat,sat)); % equation for horizontal
    slope_1(sat) = sqrt(A(1,sat)^2)*R(sat,sat)/sqrt(S(sat,sat));
    slope_2(sat) = sqrt(A(2,sat)^2)*R(sat,sat)/sqrt(S(sat,sat));
    slope_3(sat) = sqrt(A(3,sat)^2)*R(sat,sat)/sqrt(S(sat,sat));
end

% Protection level calculation
PL2 = zeros(1,3);
PL2(1) = max(slope_1)*sqrt(lambda);
PL2(2) = max(slope_2)*sqrt(lambda);
PL2(3) = max(slope_3)*sqrt(lambda);

% Solution Separation
% Simplified PL calculation
PL3 = zeros(3,1);
PL_ff = sqrt(diag(P(1:3,1:3))) * norminv(1-PHMI/3); % fault free condition

PL3 = PL3 + PL_ff;

% fault cases:
for sat_excl = 1:numsat
    satlist = setdiff(1:numsat,sat_excl);
    [~,Pss,~]=Kfilter_h(v(satlist),H(satlist,:),R(satlist,satlist),X,P);

    PL_monitor = sqrt(diag(Pss(1:3,1:3))) * norminv(1-p_sat); % fault free condition

    PL3 = PL3 + PL_monitor;
end

% Solution Separation V2

%Compute probability of no fault
% pnofault = prod(1-ones(1,numsat)*p_sat)*prod(1-p_const);
pnofault = prod(1-ones(1,numsat)*p_sat); % cannot monitor constellation
% wide fault: let's assume there is a different monitor that monitors or
% prevents this fault.

%Initialize pnotmonitored
p_not_monitored = 1-pnofault;

% determine the number of combined satellite failures that needs to be
% monitored

Nfm = 0;
u = p_sat*numsat+p_const;
for r = 1:10
    if (factorial(r)*P_THRES)^(1/r) < u & u<= (factorial(r+1)*P_THRES)^(1/(r+1))
        Nfm = r;
        break;
    end
end

% define subsets
excl_sats = nchoosek(1:numsat,Nfm);
subsets = zeros(size(excl_sats,1),numsat);
for ss = 1:size(excl_sats,1)
    for n = 1:Nfm
        subsets(ss,excl_sats(ss,n))=1;
    end
end
pap_subsets = prod((1-subsets)*(1-p_sat)+subsets).*prod(subsets*p_sat+(1-subsets));

p_not_monitored = p_not_monitored - sum(pap_subsets);

if p_not_monitored > P_THRES
    warning('P_non_monitored > P_tresh');
end

% All-in-view covariance and position
AIV_pos = X(1:3);
AIV_cov = diag(P(1:3,1:3));

% fault cases:
num_subset = numel(pap_subsets); % simplified number of subsets

PL_ss = zeros(num_subset,3);
for ss = 1:num_subset
    satlist = setdiff(1:numsat,find(subsets(ss,:)));
    [Xss,Pss,~]=Kfilter_h(v(satlist),H(satlist,:),R(satlist,satlist),X,P);

    % difference of variance
    sigma2_ss = AIV_cov - diag(Pss(1:3,1:3));

    % subset treshold
    T_ss= - norminv((1/6)*P_fa/(num_subset-1))*sqrt(abs(sigma2_ss));

    % fault-detection
    if Xss(1:3) - AIV_pos > T_ss
        warning('FDE; exclusion should be attempted')
    end

    % subset protection level
    PL_ss(ss,:) = T_ss - norminv((1/3)*PHMI / ((num_subset-1) * pap_subsets(ss))) * sqrt(diag(Pss(1:3,1:3)));
end
PL4 = max(PL_ss);

end
