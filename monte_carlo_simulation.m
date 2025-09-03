
clc; clear; close all;
N_list = 1e7; % sample number
epsilon = 0.02; % 2% instrumental error

% real value
C_Li_0_real = 5.00;
C_Na_0_real = 500.0;
C_Li_real = 3.00;
C_Na_real = 495.0;

q_Li_real = C_Li_0_real - C_Li_real;
q_Na_real = C_Na_0_real - C_Na_real;

%% 
for i = 1:length(N_list)
    N_sample = N_list(i);
    rng('default'); 
    
    % measured value
    C_Li_0_sample = normrnd(C_Li_0_real, epsilon * C_Li_0_real, [N_sample,1]);
    C_Li_sample    = normrnd(C_Li_real,    epsilon * C_Li_real,    [N_sample,1]);
    C_Na_0_sample = normrnd(C_Na_0_real, epsilon * C_Na_0_real, [N_sample,1]);
    C_Na_sample    = normrnd(C_Na_real,    epsilon * C_Na_real,    [N_sample,1]);
    
    % adsorption
    q_Li_sample = normrnd(q_Li_real, epsilon * q_Li_real, [N_sample,1]);
    q_Na_sample = normrnd(q_Na_real, epsilon * q_Na_real, [N_sample,1]);
    
    % S1
    R_Li_sample = (C_Li_0_sample - C_Li_sample) ./ C_Li_0_sample;
    R_Na_sample = (C_Na_0_sample - C_Na_sample) ./ C_Na_0_sample;
    S1 = R_Li_sample ./ R_Na_sample;
    
    % S2
    S2 = (q_Li_sample ./ C_Li_0_sample) ./ (q_Na_sample ./ C_Na_0_sample);
    
    % statistical result
    mean_S1 = mean(S1);
    std_S1 = std(S1);
    mean_S2 = mean(S2);
    std_S2 = std(S2);

    mean_R_Li=mean(R_Li_sample);
    std_R_Li=std(R_Li_sample);
    mean_R_Na=mean(R_Na_sample);
    std_R_Na=std(R_Na_sample);
    
    fprintf('N = %.0e\n', N_sample);
    fprintf('Method 1 (R-based): Mean = %.4f, Std = %.4f\n', mean_S1, std_S1);
    fprintf('Method 2 (q/C0-based): Mean = %.4f, Std = %.4f\n', mean_S2, std_S2);
    fprintf('R_Li: Mean = %.4f, Std = %.4f\n', mean_R_Li, std_R_Li);
    fprintf('R_Na: Mean = %.4f, Std = %.4f\n', mean_R_Na, std_R_Na);
    

S1_valid = S1(~isinf(S1) & ~isnan(S1) & S1 > -1e3 & S1 < 1e3); 
S2_valid = S2(~isinf(S2) & ~isnan(S2) & S2 > -1e3 & S2 < 1e3); 


figure('Name', sprintf('Monte Carlo S1 vs S2 (Histogram freq), N = %.0e', N_sample), 'Position', [100,100,800,600]);

% S1 bin:Bin k: edges(k) ≤ x < edges(k+1)

bin_width_S1 = 0.1;
edges_S1 = -200:bin_width_S1:200;

% S2 
bin_width_S2 = 0.1;
edges_S2 = -200:bin_width_S2:200;

% === S1 ===
subplot(2,1,1);
counts_S1 = histcounts(S1_valid, edges_S1);
percent_S1 = counts_S1 / length(S1_valid) * 100;
bar(edges_S1(1:end-1), percent_S1, 'FaceColor','b', 'EdgeColor','black', 'BarWidth',1);
title(sprintf('S1 (R-based) Frequency, N = %.0e', N_sample));
xlabel('S1'); ylabel('Frequency (%)');
grid on;
xline(40, 'r--', 'S=40');
xlim([-200 200]);  % 

% === S2 ===
subplot(2,1,2);
counts_S2 = histcounts(S2_valid, edges_S2);
percent_S2 = counts_S2 / length(S2_valid) * 100;
bar(edges_S2(1:end-1), percent_S2, 'FaceColor','k', 'EdgeColor','black', 'BarWidth',1);
title(sprintf('S2 (q/C0-based) Frequency, N = %.0e', N_sample));
xlabel('S2'); ylabel('Frequency (%)');
grid on;
xline(40, 'r--', 'S=40');
xlim([-200 200]);

  
    pause(0.5); % 
end
