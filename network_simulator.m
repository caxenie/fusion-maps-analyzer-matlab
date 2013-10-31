%% Fusion network for testing simple algebraic and temporal relationships
% cleanup
close all;
clear all;
clc;

%% INITIALIZATION
%% Simulation parameters
% default values for learning rate
ETA  = 0.002;
ETAH = 10*ETA;
% add noise in the sensor signal
noise_on = 0;
% iterations to present a data point to the net
show_step          = 1000;
% prepare sensor data to fed to the net
run_steps          = 20*show_step;
% runtime parameters
run_iters          = run_steps - 1;
run_iters_extended = run_iters + 1;
% number of maps
maps               = 4;
% complex / simple learning rate flag
learning_update_type = 'simple'; % {simple, complex}
% check learning rate adaptation rule
switch learning_update_type
    case 'simple'
        % error type: simple difference
        % {fixed, adaptive, decay, divisive}
        type = 'fixed';
    case 'complex'
        % error type: squared difference
        % {grad-history, grad-history-average}
        type = 'grad-history';
end
% id of the current update rule
rand_update_rule  = 1;
% possible update rules - depends on net topology
update_rules      = 9;
% individual id for each rule
rules_ids         = 1:update_rules;
% convergence steps counter
convergence_steps = 1;
% simulation points (sensor data is presented for half simulation)
sim_points        = run_iters_extended - 1;
%% Network elements setup
% sensor values init - each sensor has it's own model
m1_sensor = rand+zeros(1, run_iters);
m2_sensor = rand+zeros(1, run_iters);
m3_sensor = rand+zeros(1, run_iters);
m4_sensor = rand+zeros(1, run_iters);

%% input pattern of the sensors (TODO make it look nicer)
increment = 0.0001;
switch1_s1 = 3000;
switch2_s1 = 5000;
% switch1_s2 = 2000;
% switch2_s2 = 6000;
switch1_s3 = 10000;
switch2_s3 = 12000;
switch3_s3 = 16000;
% switch1_s4 = 16000;
% switch2_s4 = 17000;
% switch3_s4 = 18000;
% sensor 1
for i=2:run_iters
    if(i<=switch1_s1)
        m1_sensor(i) = m1_sensor(i-1) + increment;
    end
    if(i>switch1_s1 && i<= switch2_s1)
        m1_sensor(i) = m1_sensor(i-1);
    end
end
% % sensor 2
% for i=2:run_iters
%     if(i<=switch1_s2)
%         m2_sensor(i) = m2_sensor(i-1) + increment;
%     end
%     if(i>switch1_s2 && i<= switch2_s2)
%         m2_sensor(i) = m2_sensor(i-1);
%     end
% end
% sensor 3
for i=2:run_iters
    if(i>=switch1_s3 && i<= switch2_s3)
        m3_sensor(i) = m3_sensor(i-1) + increment;
    end
    if(i>switch2_s3 && i<= switch3_s3)
        m3_sensor(i) = m3_sensor(i-1);
    end
end
% sensor 4
% for i=2:run_iters
%     if(i>=switch1_s4 && i<=switch2_s4)
%         m4_sensor(i) = m4_sensor(i-1) + increment;
%     end
%     if(i>switch2_s4 && i<= switch3_s4)
%         m4_sensor(i) = m4_sensor(i-1);
%     end
% end

% add noise over the sensor signal
if(noise_on==1)
    sigma = 0.0025; % noise standard deviation
    m1_sensor = m1_sensor + sigma*randn(size(m1_sensor));
    m2_sensor = m2_sensor + sigma*randn(size(m2_sensor));
    m3_sensor = m3_sensor + sigma*randn(size(m3_sensor));
    m4_sensor = m4_sensor + sigma*randn(size(m4_sensor));
end

% network iterator
net_iter = 1;

% init maps and indices
m1 = rand; m2 = rand; m3 = rand; m4 = rand;

% maps ids
m1_id = 1; m2_id = 2; m3_id = 3; m4_id = 4;

% number of connections of the maps (sensor + relations to other maps)
m1_links = 1:2; m2_links = 1:3;
m3_links = 1:2; m4_links = 1:2;

% init maps errors for learning rate adaptation
em1 = zeros(1, length(m1_links)); em2 = zeros(1, length(m2_links));
em3 = zeros(1, length(m3_links)); em4 = zeros(1, length(m4_links));

% error gradients w.r.t. each map initialization (for learning rate adaptation)
dem1 = zeros(1, length(m1_links)); dem2 = zeros(1, length(m2_links));
dem3 = zeros(1, length(m3_links)); dem4 = zeros(1, length(m4_links));

dem1_old = zeros(1, length(m1_links)); dem2_old = zeros(1, length(m2_links));
dem3_old = zeros(1, length(m3_links)); dem4_old = zeros(1, length(m4_links));

% global error gradient history all maps errors
grad_e = zeros(run_iters_extended, maps);

% delta-bar-delta learning rate adaptation
grad_bar1 = zeros(1, length(m1_links)); grad_bar2 = zeros(1, length(m2_links));
grad_bar3 = zeros(1, length(m3_links)); grad_bar4 = zeros(1, length(m4_links));

grad_bar1_old = zeros(1, length(m1_links)); grad_bar2_old = zeros(1, length(m2_links));
grad_bar3_old = zeros(1, length(m3_links)); grad_bar4_old = zeros(1, length(m4_links));

% learning rates setup
maps_nr     = maps;
sensors_nr  = maps_nr;
error_nr    = update_rules;
lrates_nr   = update_rules;
net_data    = zeros(sim_points, maps_nr+error_nr+lrates_nr);

% init learning rates and bounds
etam1 = ETA*ones(sim_points, length(m1_links));
etam2 = ETA*ones(sim_points, length(m2_links));
etam3 = ETA*ones(sim_points, length(m3_links));
etam4 = ETA*ones(sim_points, length(m4_links));

%% NETWORK DYNAMICS
%
% encoded relationships
%
%   m2 = 3*m1
%   m3 = m2*m4
%
%   noticeable changes from baseline in m1 and m3
%

while(1)
    % run the net to prove convergence after sensor decoupling
    if (convergence_steps == run_iters_extended)
        break;
    end;
       
    % shuffle maps ids for update
    idx = update_rules;
    while idx>1
        jdx = randi(update_rules, 1);
        tmp = rules_ids(jdx);
        rules_ids(jdx) = rules_ids(idx);
        rules_ids(idx) = tmp;
        idx = idx - 1;
    end
    
    %% go through current convergence step update sequence
    if(strcmp(learning_update_type,'simple')==1)
        for i = 1:update_rules
            % pick a rule to update
            rand_update_rule = rules_ids(i);
            % update the corresponding node using simple error
            switch rand_update_rule
                % Circular permutations of the error signal
                %% m1 updates
                case 1
                    m1 = (1-etam1(convergence_steps, m1_links(1)))*m1 + etam1(convergence_steps, m1_links(1))*...
                        m1_sensor(net_iter);
                case 2
                    m1 = (1-etam1(convergence_steps, m1_links(2)))*m1 + etam1(convergence_steps, m1_links(2)) * ...
                        1/3*m2;
                    %% m2 updates
                case 3
                    m2 = (1-etam2(convergence_steps, m2_links(1)))*m2 + etam2(convergence_steps, m2_links(1))*...
                        m2_sensor(net_iter);
                case 4
                    m2 = (1-etam2(convergence_steps, m2_links(2)))*m2 + etam2(convergence_steps, m2_links(2)) * ...
                        3*m1;
                case 5
                    m2 = (1-etam2(convergence_steps, m2_links(3)))*m2 + etam2(convergence_steps, m2_links(3)) * ...
                        m3/m4;
                    %% m3 updates
                case 6
                    m3 = (1-etam3(convergence_steps, m3_links(1)))*m3 + etam3(convergence_steps, m3_links(1))*...
                        m3_sensor(net_iter);
                case 7
                    m3 = (1-etam3(convergence_steps, m3_links(2)))*m3 + etam3(convergence_steps, m3_links(2)) * ...
                        m2*m4;
                    %% m4 updates
                case 8
                    m4 = (1-etam4(convergence_steps, m4_links(1)))*m4 + etam4(convergence_steps, m4_links(1))*...
                        m4_sensor(net_iter);
                case 9
                    m4 = (1-etam4(convergence_steps, m4_links(2)))*m4 + etam4(convergence_steps, m4_links(2)) * ...
                        m3/m2;
            end
        end
    else
        for i = 1:update_rules
            % pick a rule to update
            rand_update_rule = rules_ids(i);
            % the complex update rules (squared error)
            % update the corresponding node
            switch rand_update_rule
                %% m1 updates
                case 1
                    m1 = (1-2*etam1(convergence_steps, m1_links(1)))*m1 + 2*etam1(convergence_steps, m1_links(1))*...
                        m1_sensor(net_iter);
                case 2
                    m1 = (1-18*etam1(convergence_steps, m1_links(2)))*m1 + 6*etam1(convergence_steps, m1_links(2)) * ...
                        m2;
                    %% m2 updates
                case 3
                    m2 = (1-2*etam2(convergence_steps, m2_links(1)))*m2 + 2*etam2(convergence_steps, m2_links(1))*...
                        m2_sensor(net_iter);
                case 4
                    m2 = (1-2*etam2(convergence_steps, m2_links(2)))*m2 + 6*etam2(convergence_steps, m2_links(2)) * ...
                        m1;
                case 5
                    m2 = (1-2*m4^2*etam2(convergence_steps, m2_links(3)))*m2 + 2*etam2(convergence_steps, m2_links(3)) * ...
                        m3*m4;
                    %% m3 updates
                case 6
                    m3 = (1-2*etam3(convergence_steps, m3_links(1)))*m3 + 2*etam3(convergence_steps, m3_links(1))*...
                        m3_sensor(net_iter);
                case 7
                    m3 = (1-2*etam3(convergence_steps, m3_links(2)))*m3 + 2*etam3(convergence_steps, m3_links(2)) * ...
                        m2*m4;
                    %% m4 updates
                case 8
                    m4 = (1-2*etam4(convergence_steps, m4_links(1)))*m4 + 2*etam4(convergence_steps, m4_links(1))*...
                        m4_sensor(net_iter);
                case 9
                    m4 = (1-2*m2^2*etam4(convergence_steps, m4_links(2)))*m4 + 2*etam4(convergence_steps, m4_links(2)) * ...
                        m3*m2;
            end
        end
        
    end
    
    %% LEARNING RATES ADAPTATION
    
    switch(learning_update_type)
        case 'simple'
            % --------errors--------
            % error computation for learning rate adaptation
            % m1
            em1(1) = m1 - m1_sensor(net_iter);
            em1(2) = m1 - m2/3;
            % m2
            em2(1) = m2 - m2_sensor(net_iter);
            em2(2) = m2 - 3*m1;
            em2(3) = m2 - m3/m4;
            % m3
            em3(1) = m3 - m3_sensor(net_iter);
            em3(2) = m3 - m2*m4;
            % m4
            em4(1) = m4 - m4_sensor(net_iter);
            em4(2) = m4 - m3/m2;
            
            % get only the amplitude of the error
            em1(1) = abs(em1(1));
            em1(2) = abs(em1(2));
            em2(1) = abs(em2(1));
            em2(2) = abs(em2(2));
            em2(3) = abs(em2(3));
            em3(1) = abs(em3(1));
            em3(2) = abs(em3(2));
            em4(1) = abs(em4(1));
            em4(2) = abs(em4(2));
            
            % simple update rules for learning rate adaptation (divisive/decay)
            for k=1:length(m1_links)
                etam1(convergence_steps, m1_links(k)) = update_learning_rate(etam1(convergence_steps, m1_links(k)), em1, m1_links(k), ETA, type);
                etam1(convergence_steps, m1_links(k)) = clamp(etam1(convergence_steps, m1_links(k)), ETAH);
            end
            
            for k=1:length(m2_links)
                etam2(convergence_steps, m2_links(k)) = update_learning_rate(etam2(convergence_steps, m2_links(k)), em2, m2_links(k), ETA, type);
                etam2(convergence_steps, m2_links(k)) = clamp(etam2(convergence_steps, m2_links(k)), ETAH);
            end
            for k=1:length(m3_links)
                etam3(convergence_steps, m3_links(k)) = update_learning_rate(etam3(convergence_steps, m3_links(k)), em3, m3_links(k), ETA, type);
                etam3(convergence_steps, m3_links(k)) = clamp(etam3(convergence_steps, m3_links(k)), ETAH);
            end
            for k=1:length(m4_links)
                etam4(convergence_steps, m4_links(k)) = update_learning_rate(etam4(convergence_steps, m4_links(k)), em4, m4_links(k), ETA, type);
                etam4(convergence_steps, m4_links(k)) = clamp(etam4(convergence_steps, m4_links(k)), ETAH);
            end
            
        case 'complex'
            % -----------errors------------
            % error computation for learning rate adaptation (squared err)
            % m1
            em1(1) = (m1 - m1_sensor(net_iter))^2;
            em1(2) = (m2 - 3*m1)^2;
            % m2
            em2(1) = (m2 - m2_sensor(net_iter))^2;
            em2(2) = (m2 - 3*m1)^2;
            em2(3) = (m3 - m2*m4)^2;
            % m3
            em3(1) = (m3 - m3_sensor(net_iter))^2;
            em3(2) = (m3 -m2*m4)^2;
            % m4
            em4(1) = (m4 - m4_sensor(net_iter))^2;
            em4(2) = (m3 - m2*m4)^2;
            
            % ---------gradients----------
            % gradient of errors for m1
            dem1(m1_links(1)) = 2*(m1-m1_sensor(net_iter))*(-1);
            dem1(m1_links(2)) = 2*(m2-3*m1)*(-3);
            
            % gradient of errors for m2
            dem2(m2_links(1)) = 2*(m2-m2_sensor(net_iter))*(-1);
            dem2(m2_links(2)) = 2*(m2-3*m1);
            dem2(m2_links(3)) = 2*(m3 - m2*m4)*(-m4);
            
            % gradient of errors for m3
            dem3(m3_links(1)) = 2*(m3-m3_sensor(net_iter))*(-1);
            dem3(m3_links(2)) = 2*(m3-m2*m4);
            
            % gradient of errors for m4
            dem4(m4_links(1)) = 2*(m4-m4_sensor(net_iter))*(-1);
            dem4(m4_links(2)) = 2*(m3-m2*m4)*(-m2);
            
            % params init for grad-history
            u       = 1.6;            % u > 1
            d       = 1/u;            % d < 1
            l_min   = ETA;
            l_max   = ETAH;
            
            % params for grad-history-average
            beta    = 0.004;         % 0 < beta < 1
            k       = 0.004;
            gamma   = 0.8;
            
            % complex update rules for the learning rate
            for k=1:length(m1_links)
                % bar grad
                grad_bar1(k) = (1-beta)*dem1(k) + beta*grad_bar1_old(k);
                % learning rates5
                etam1(convergence_steps, m1_links(k)) = update_learning_rate_complex(etam1(convergence_steps, m1_links(k)),  grad_bar1_old(k),  dem1(k), dem1_old(k), u, d, l_min, l_max, k, gamma, type);
                etam1(convergence_steps, m1_links(k)) = clamp(etam1(convergence_steps, m1_links(k)), ETAH);
            end
            for k=1:length(m2_links)
                % bar grad
                grad_bar2(k) = (1-beta)*dem2(k) + beta*grad_bar2_old(k);
                % learning rates
                etam2(convergence_steps, m2_links(k)) = update_learning_rate_complex(etam2(convergence_steps, m2_links(k)),  grad_bar2_old(k),  dem2(k), dem2_old(k), u, d, l_min, l_max, k, gamma, type);
                etam2(convergence_steps, m2_links(k)) = clamp(etam2(convergence_steps, m2_links(k)), ETAH);
            end
            
            for k=1:length(m3_links)
                % bar grad
                grad_bar3(k) = (1-beta)*dem3(k) + beta*grad_bar3_old(k);
                % learning rates
                etam3(convergence_steps, m3_links(k)) = update_learning_rate_complex(etam3(convergence_steps, m3_links(k)),  grad_bar3_old(k),  dem3(k), dem3_old(k), u, d, l_min, l_max, k, gamma, type);
                etam3(convergence_steps, m3_links(k)) = clamp(etam3(convergence_steps, m3_links(k)), ETAH);
            end
            
            for k=1:length(m4_links)
                % bar grad
                grad_bar4(k) = (1-beta)*dem4(k) + beta*grad_bar4_old(k);
                % learning rates
                etam4(convergence_steps, m4_links(k)) = update_learning_rate_complex(etam4(convergence_steps, m4_links(k)),  grad_bar4_old(k),  dem4(k), dem4_old(k), u, d, l_min, l_max, k, gamma, type);
                etam4(convergence_steps, m4_links(k)) = clamp(etam4(convergence_steps, m4_links(k)), ETAH);
            end
            
            % update indices
            dem1_old = dem1; dem2_old = dem2; dem3_old = dem3; dem4_old = dem4;
            grad_bar1_old = grad_bar1;
            grad_bar2_old = grad_bar2;
            grad_bar3_old = grad_bar3;
            grad_bar4_old = grad_bar4;
            
            % error gradient history for analysis
            grad_e(convergence_steps, 1) = dem1(1);
            grad_e(convergence_steps, 2) = dem1(2);
            grad_e(convergence_steps, 3) = dem2(1);
            grad_e(convergence_steps, 4) = dem2(2);
            grad_e(convergence_steps, 5) = dem2(3);
            grad_e(convergence_steps, 6) = dem3(1);
            grad_e(convergence_steps, 7) = dem3(2);
            grad_e(convergence_steps, 8) = dem4(1);
            grad_e(convergence_steps, 9) = dem4(2);
    end
    
    %% WRITE DATA TO STRUCT
    % maps
    net_data(convergence_steps, 1) = m1;
    net_data(convergence_steps, 2) = m2;
    net_data(convergence_steps, 3) = m3;
    net_data(convergence_steps, 4) = m4;
    % error
    net_data(convergence_steps, 7) = em1(1);
    net_data(convergence_steps, 8) = em1(2);
    
    net_data(convergence_steps, 9)  = em2(1);
    net_data(convergence_steps, 10) = em2(2);
    net_data(convergence_steps, 11) = em2(3);
    
    net_data(convergence_steps, 12) = em3(1);
    net_data(convergence_steps, 13) = em3(2);
    
    net_data(convergence_steps, 14) = em4(1);
    net_data(convergence_steps, 15) = em4(2);
    
    % sample idx
    net_data(convergence_steps, 21) = convergence_steps ;
    
    % add the learning rates in the struct
    net_data(convergence_steps, 22) = etam1(convergence_steps, m1_links(1));
    net_data(convergence_steps, 23) = etam1(convergence_steps, m1_links(2));
    
    net_data(convergence_steps, 24) = etam2(convergence_steps, m2_links(1));
    net_data(convergence_steps, 25) = etam2(convergence_steps, m2_links(2));
    net_data(convergence_steps, 26) = etam2(convergence_steps, m2_links(3));
    
    net_data(convergence_steps, 27) = etam3(convergence_steps, m3_links(1));
    net_data(convergence_steps, 28) = etam3(convergence_steps, m3_links(2));
    
    net_data(convergence_steps, 29) = etam4(convergence_steps, m4_links(1));
    net_data(convergence_steps, 30) = etam4(convergence_steps, m4_links(2));
    
    %% update indices
    net_iter = net_iter + 1;
    convergence_steps = convergence_steps + 1;
end

% fill in the net simulation data into the visualization struct
fusion_analyzer_data = net_data;

%% VISUALIZATION
figure(1);
% map 1 and sensor 1
hd(1) = subplot(2, 2, 1);
plot(net_data(:, 1), '.r'); hold on;
plot(m1_sensor, '.k');
grid on; legend('M1 map', 'Sensor 1');
% map 2 and sensor 2
hd(2) = subplot(2, 2, 2);
plot(net_data(:, 2), '.b'); hold on;
plot(m2_sensor, '.k');
grid on; legend('M2 map','Sensor 2');
% map 3 and sensor 3
hd(5) = subplot(2, 2, 3);
plot(net_data(:, 3), '.g'); hold on;
plot(m3_sensor, '.k');
grid on; legend('M3 map ', 'Sensor 3');
% map 4 and sensor 4
hd(6) = subplot(2, 2, 4);
plot(net_data(:, 4), '.c'); hold on;
plot(m4_sensor, '.k');
grid on; legend('M4 map', 'Sensor 4');

% link axes
linkaxes(hd, 'x');
% set figure props
set(gcf,'color','w');

% -------------- Errors signals --------------
figure(2);
%-------------------M1-------------------
he(1) = subplot(2,2,1);
plot(net_data(:,7),'.k');hold on;
plot(net_data(:,8), '.m');
title('M1 errors');
legend('Err w.r.t S1', 'Err w.r.t R1');
grid on;
%-------------------M2-------------------
he(2) = subplot(2,2,2);
plot(net_data(:,9), '.k'); hold on;
plot(net_data(:,10), '.m'); hold on;
plot(net_data(:,11), '.y');
title('M2 errors');
legend('Err w.r.t S2', 'Err w.r.t R1', 'Err w.r.t R2');
grid on;
%-------------------M3-------------------
he(3) = subplot(2,2,3);
plot(net_data(:,12), '.k'); hold on;
plot(net_data(:, 13), '.y');
title('M3 errors');
legend('Err w.r.t S3', 'Err w.r.t R2');
grid on;
%-------------------M4-------------------
he(4) = subplot(2,2,4);
plot(net_data(:,14), '.k'); hold on;
plot(net_data(:,15), '.y');
title('M4 errors');
legend('Err w.r.t S4', 'Err w.r.t R2');
grid on;

% set figure props
set(gcf,'color','w');

% plot learning rates on a per map basis
figure(3);
% ------------------m1-----------------
heta1 = subplot(2,2,1);
plot(fusion_analyzer_data(:, 22), '.k'); hold on;
plot(fusion_analyzer_data(:, 23), '.m');
grid on;
legend('Lrate w.r.t S1','Lrate w.r.t. R1');
title('Learning rate analysis');
% ------------------m2-----------------
heta2 = subplot(2,2,2);
plot(fusion_analyzer_data(:, 24), '.k'); hold on;
plot(fusion_analyzer_data(:, 25), '.m'); hold on;
plot(fusion_analyzer_data(:, 26), '.y');
grid on; legend('Lrate w.r.t S2','Lrate w.r.t. R1', 'Lrate w.r.t. R2');
% ------------------m3-----------------
heta3 = subplot(2,2,3);
plot(fusion_analyzer_data(:, 27), '.k'); hold on;
plot(fusion_analyzer_data(:, 28), '.y');
grid on; legend('Lrate w.r.t S3','Lrate w.r.t. R2');
% ------------------m4-----------------
heta4 = subplot(2,2,4);
plot(fusion_analyzer_data(:, 29), '.k'); hold on;
plot(fusion_analyzer_data(:, 30), '.y');
grid on; legend('Lrate w.r.t S4','Lrate w.r.t. R2' );

heta = [heta1 heta2 heta3 heta4];
% link axes for analysis
linkaxes(heta, 'x');

% set figure props
set(gcf,'color','w');

%% plot the error gradients
if(strcmp(learning_update_type,'complex')==1)
    figure;
    set(gcf,'color','w');
    % -------m1--------
    subplot(2,3,1);
    plot(grad_e(:,1),'.k'); hold on;
    plot(grad_e(:,2),'.m');
    title('m1 errors derivatives'); grid on;
    legend('dE w.r.t. sensor','dE w.r.t. relation 1');
    % -------m2--------
    subplot(2,3,2);
    plot(grad_e(:,3),'.k'); hold on;
    plot(grad_e(:,4),'.m');hold on;
    plot(grad_e(:,5), '.y');
    title('m2 errors derivatives'); grid on;
    legend('dE w.r.t. sensor','dE w.r.t. relation 1','dE w.r.t. relation 2');
    % -------m3--------
    subplot(2,3,3);
    plot(grad_e(:,6),'.k'); hold on;
    plot(grad_e(:,7),'.y');
    title('m3 errors derivatives'); grid on;
    legend('dE w.r.t. sensor','dE w.r.t. relation 2');
    % -------m4--------
    subplot(2,3,4);
    plot(grad_e(:,8),'.k'); hold on;
    plot(grad_e(:,9),'.y');
    title('m4 errors derivatives'); grid on;
    legend('dE w.r.t. sensor','dE w.r.t. relation 2');
    % set figure props
    set(gcf,'color','w');
    
end


% plot the error profiles 
figure
% ----- m1 ------
subplot(5, 2, 1);
plot(net_data(:,8)./net_data(:,7), '.r');
title('Err w.r.t R1/Err w.r.t S');
grid on;
subplot(5, 2, 3);
plot(net_data(:,7)./net_data(:,8), '.r');
title('Err w.r.t S/Err w.r.t R1');
grid on;
% ----- m2 ------
subplot(5, 2, 2);
plot((net_data(:,10) + net_data(:,11))./net_data(:,9), '.b');
title('(Err w.r.t R1+Err w.r.t R2)/Err w.r.t S');
grid on;
subplot(5, 2, 4);
plot((net_data(:,9) + net_data(:,11))./net_data(:,10), '.b');
title('(Err w.r.t S+Err w.r.t R1)/Err w.r.t R1');
grid on;
subplot(5, 2, 6);
plot((net_data(:,9) + net_data(:,10))./net_data(:,11), '.b');
title('(Err w.r.t S+Err w.r.t R2)/Err w.r.t R2');
grid on;
% ----- m3 ------
subplot(5, 2, 7);
plot(net_data(:,13)./net_data(:,12), '.g');
title('Err w.r.t R2/Err w.r.t S');
grid on;
subplot(5, 2, 9);
plot(net_data(:,12)./net_data(:,13), '.g');
title('Err w.r.t S/Err w.r.t R2');
grid on;
% ----- m4 ------
subplot(5, 2, 8);
plot(net_data(:,15)./net_data(:,14), '.y');
title('Err w.r.t R2/Err w.r.t S');
grid on;
subplot(5, 2, 10);
plot(net_data(:,14)./net_data(:,15), '.y');
title('Err w.r.t S/Err w.r.t R2');
grid on;
% set figure props
set(gcf,'color','w');