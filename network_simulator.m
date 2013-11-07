%% Fusion network for testing simple algebraic and temporal relationships
% cleanup
close all;
clear all;
clc;

%% INITIALIZATION
%% Simulation parameters
% default values for confidence factor
ETA  = 0.002;
ETAH = 10*ETA;
% add noise in the sensor signal
noise_on = 0;
% prepare sensor data to fed to the net
run_steps          = 40000;
% runtime parameters
run_iters          = run_steps - 1;
run_iters_extended = run_iters + 1;
% number of maps
maps               = 4;
% complex / simple confidence factor selector flag
confidence_factor_type = 'simple'; % {simple, complex}
% check confidence factor adaptation rule
switch confidence_factor_type
    case 'simple'
        % error type: simple difference
        % {fixed, adaptive, decay, divisive}
        type = 'decay';
    case 'complex'
        % error type: squared difference
        % {grad-history, grad-history-average}
        type = 'grad-history-average';
end
% id of the current update rule
rand_update_rule  = 1;
% possible update rules - depends on net topology
update_rules      = 9;
% individual id for each rule
rules_ids         = 1:update_rules;
% simulation points (sensor data is presented for half simulation)
sim_points        = run_iters_extended - 1;
%% Network elements setup
% sensor values init (baseline for each sensor so that net is relaxed)
m1_sensor = 0.5  + zeros(1, run_iters);
m2_sensor = 1.5  + zeros(1, run_iters);
m3_sensor = 0.3  + zeros(1, run_iters);
m4_sensor = 0.2  + zeros(1, run_iters);

% convergence steps counter
convergence_steps = 1;

%% input pattern of the sensors (TODO make it look nicer)
increment = 0.0001;
switch1_s1 = 10000;
switch2_s1 = 12000;
switch3_s1 = 20000;
% switch1_s2 = 2000;
% switch2_s2 = 6000;
switch1_s3 = 26000;
switch2_s3 = 29000;
switch3_s3 = 35000;
% switch1_s4 = 16000;
% switch2_s4 = 17000;
% switch3_s4 = 18000;
% sensor 1
for i=2:run_iters
    if(i>=switch1_s1 && i<=switch2_s1)
        m1_sensor(i) = m1_sensor(i-1) + increment;
    end
    if(i>switch2_s1 && i<= switch3_s1)
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
        m3_sensor(i) = m3_sensor(i-1) - increment;
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

% init maps errors for confidence factor adaptation
em1 = zeros(1, length(m1_links)); em2 = zeros(1, length(m2_links));
em3 = zeros(1, length(m3_links)); em4 = zeros(1, length(m4_links));

% error gradients w.r.t. each map initialization (for confidence factor adaptation)
dem1 = zeros(1, length(m1_links)); dem2 = zeros(1, length(m2_links));
dem3 = zeros(1, length(m3_links)); dem4 = zeros(1, length(m4_links));

dem1_old = zeros(1, length(m1_links)); dem2_old = zeros(1, length(m2_links));
dem3_old = zeros(1, length(m3_links)); dem4_old = zeros(1, length(m4_links));

% global error gradient history all maps errors
grad_e = zeros(run_iters_extended, maps);

% delta-bar-delta confidence factor adaptation
grad_bar1 = zeros(1, length(m1_links)); grad_bar2 = zeros(1, length(m2_links));
grad_bar3 = zeros(1, length(m3_links)); grad_bar4 = zeros(1, length(m4_links));

grad_bar1_old = zeros(1, length(m1_links)); grad_bar2_old = zeros(1, length(m2_links));
grad_bar3_old = zeros(1, length(m3_links)); grad_bar4_old = zeros(1, length(m4_links));

% confidence factors setup
maps_nr     = maps;
sensors_nr  = maps_nr;
error_nr    = update_rules;
lrates_nr   = update_rules;
net_data    = zeros(sim_points, maps_nr+error_nr+lrates_nr);

% init confidence factors and bounds
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
%   noticeable changes from sensor baseline in m1 and m3
%

while(1)
    % end of simulation 
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
       
    %% COMPUTE UPDATE SEQUENCE IN THE CURRENT CONVERGENCE STEP
    if(strcmp(confidence_factor_type,'simple')==1)
        for i = 1:update_rules
            % pick a rule to update
            rand_update_rule = rules_ids(i);
            
                % update the corresponding node using simple error
                switch rand_update_rule
                    % Circular permutations of the error signal
                    %% m1 updates
                    case 1
                        m1 = m1 - etam1(convergence_steps, m1_links(1))*em1(1);
                    case 2
                        m1 = m1 - etam1(convergence_steps, m1_links(2))*em1(2);
                        %% m2 updates
                    case 3
                        m2 = m2 - etam2(convergence_steps, m2_links(1))*em2(1);
                    case 4
                        m2 = m2 - etam2(convergence_steps, m2_links(2))*em2(2);
                    case 5
                        m2 = m2 - etam2(convergence_steps, m2_links(3))*em2(3);
                        %% m3 updates
                    case 6
                        m3 = m3 - etam3(convergence_steps, m3_links(1))*em3(1);
                    case 7
                        m3 = m3 - etam3(convergence_steps, m3_links(2))*em3(2);
                        %% m4 updates
                    case 8
                        m4 = m4 - etam4(convergence_steps, m4_links(1))*em4(1);
                    case 9
                        m4 = m4 - etam4(convergence_steps, m4_links(2))*em4(2);
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
    
    %% COMPUTE LOCAL MISMATCH USING ONLY LOCAL INFO - CONFIDENCE ADAPTATION
    
    switch(confidence_factor_type)
        case 'simple'
            % --------errors--------
            % error computation for confidence factor adaptation
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
           
            % check if the update rules are using the divisive factor in
            % their expression
            if(strcmp(type, 'divisive')==1 || strcmp(type, 'decay')==1)
                % simple update rules for confidence factor adaptation
                for k=1:length(m1_links)
                    etam1(convergence_steps, m1_links(k)) = update_confidence_factor(etam1(convergence_steps, m1_links(k)), em1, m1_links(k), ETA, type);
                    etam1(convergence_steps, m1_links(k)) = clamp(etam1(convergence_steps, m1_links(k)), ETAH);
                end
                
                for k=1:length(m2_links)
                    etam2(convergence_steps, m2_links(k)) = update_confidence_factor(etam2(convergence_steps, m2_links(k)), em2, m2_links(k), ETA, type);
                    etam2(convergence_steps, m2_links(k)) = clamp(etam2(convergence_steps, m2_links(k)), ETAH);
                end
                for k=1:length(m3_links)
                    etam3(convergence_steps, m3_links(k)) = update_confidence_factor(etam3(convergence_steps, m3_links(k)), em3, m3_links(k), ETA, type);
                    etam3(convergence_steps, m3_links(k)) = clamp(etam3(convergence_steps, m3_links(k)), ETAH);
                end
                for k=1:length(m4_links)
                    etam4(convergence_steps, m4_links(k)) = update_confidence_factor(etam4(convergence_steps, m4_links(k)), em4, m4_links(k), ETA, type);
                    etam4(convergence_steps, m4_links(k)) = clamp(etam4(convergence_steps, m4_links(k)), ETAH);
                end
                
                
            else
                
                % simple update rules for confidence factor adaptation
                for k=1:length(m1_links)
                    etam1(convergence_steps+1, m1_links(k)) = update_confidence_factor(etam1(convergence_steps, m1_links(k)), em1, m1_links(k), ETA, type);
                    etam1(convergence_steps+1, m1_links(k)) = clamp(etam1(convergence_steps+1, m1_links(k)), ETAH);
                end
                
                for k=1:length(m2_links)
                    etam2(convergence_steps+1, m2_links(k)) = update_confidence_factor(etam2(convergence_steps, m2_links(k)), em2, m2_links(k), ETA, type);
                    etam2(convergence_steps+1, m2_links(k)) = clamp(etam2(convergence_steps+1, m2_links(k)), ETAH);
                end
                for k=1:length(m3_links)
                    etam3(convergence_steps+1, m3_links(k)) = update_confidence_factor(etam3(convergence_steps, m3_links(k)), em3, m3_links(k), ETA, type);
                    etam3(convergence_steps+1, m3_links(k)) = clamp(etam3(convergence_steps+1, m3_links(k)), ETAH);
                end
                for k=1:length(m4_links)
                    etam4(convergence_steps+1, m4_links(k)) = update_confidence_factor(etam4(convergence_steps, m4_links(k)), em4, m4_links(k), ETA, type);
                    etam4(convergence_steps+1, m4_links(k)) = clamp(etam4(convergence_steps+1, m4_links(k)), ETAH);
                end
                
            end
            
        case 'complex'
            % -----------errors------------
            % error computation for confidence factor adaptation (squared err)
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
            
                
             % complex update rules for the confidence factor
            for k=1:length(m1_links)
                % bar grad
                grad_bar1(k) = (1-beta)*dem1(k) + beta*grad_bar1_old(k);
                % confidence factors5
                etam1(convergence_steps+1, m1_links(k)) = update_confidence_factor_complex(etam1(convergence_steps, m1_links(k)),  grad_bar1_old(k),  dem1(k), dem1_old(k), u, d, l_min, l_max, k, gamma, type);
                etam1(convergence_steps+1, m1_links(k)) = clamp(etam1(convergence_steps+1, m1_links(k)), ETAH);
            end
            for k=1:length(m2_links)
                % bar grad
                grad_bar2(k) = (1-beta)*dem2(k) + beta*grad_bar2_old(k);
                % confidence factors
                etam2(convergence_steps+1, m2_links(k)) = update_confidence_factor_complex(etam2(convergence_steps, m2_links(k)),  grad_bar2_old(k),  dem2(k), dem2_old(k), u, d, l_min, l_max, k, gamma, type);
                etam2(convergence_steps+1, m2_links(k)) = clamp(etam2(convergence_steps+1, m2_links(k)), ETAH);
            end
            
            for k=1:length(m3_links)
                % bar grad
                grad_bar3(k) = (1-beta)*dem3(k) + beta*grad_bar3_old(k);
                % confidence factors
                etam3(convergence_steps+1, m3_links(k)) = update_confidence_factor_complex(etam3(convergence_steps, m3_links(k)),  grad_bar3_old(k),  dem3(k), dem3_old(k), u, d, l_min, l_max, k, gamma, type);
                etam3(convergence_steps+1, m3_links(k)) = clamp(etam3(convergence_steps, m3_links(k)), ETAH);
            end
            
            for k=1:length(m4_links)
                % bar grad
                grad_bar4(k) = (1-beta)*dem4(k) + beta*grad_bar4_old(k);
                % confidence factors
                etam4(convergence_steps+1, m4_links(k)) = update_confidence_factor_complex(etam4(convergence_steps, m4_links(k)),  grad_bar4_old(k),  dem4(k), dem4_old(k), u, d, l_min, l_max, k, gamma, type);
                etam4(convergence_steps+1, m4_links(k)) = clamp(etam4(convergence_steps+1, m4_links(k)), ETAH);
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
        
        % add the confidence factors in the struct
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

% save data in file for comparative visualization 
dump_file_name = strcat('network_dataset_',type);
save(dump_file_name, 'net_data');

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
% plot only the error magnitudes
net_data(:,7) = abs(net_data(:,7));
net_data(:,8) = abs(net_data(:,8));
net_data(:,9) = abs(net_data(:,9));
net_data(:,10) = abs(net_data(:,10));
net_data(:,11) = abs(net_data(:,11));
net_data(:,12) = abs(net_data(:,12));
net_data(:,13) = abs(net_data(:,13));
net_data(:,14) = abs(net_data(:,14));
net_data(:,15) = abs(net_data(:,15));
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

% link axes 
linkaxes(he, 'x');

% set figure props
set(gcf,'color','w');

% for adaptive is good to see the switches in the values as well
if(strcmp(type,'adaptive')==0)
    % plot confidence factors on a per map basis
    figure(3);
    % ------------------m1-----------------
    heta1 = subplot(2,2,1);
    plot(fusion_analyzer_data(:, 22), '.k'); hold on;
    plot(fusion_analyzer_data(:, 23), '.m');
    grid on;
    legend('Lrate w.r.t S1','Lrate w.r.t. R1');
    title('confidence factor analysis');
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
else
    % plot confidence factors on a per map basis
    figure(3);
    % ------------------m1-----------------
    heta1 = subplot(2,2,1);
    h(1) = plot(fusion_analyzer_data(:, 22), '-k'); hold on;
    h(2) = plot(fusion_analyzer_data(:, 23), '-m');
    set(h,'LineWidth',4);
    grid on;
    legend('Lrate w.r.t S1','Lrate w.r.t. R1');
    title('confidence factor analysis');
    % ------------------m2-----------------
    heta2 = subplot(2,2,2);
    h(1) = plot(fusion_analyzer_data(:, 24), '-k'); hold on;
    h(2) = plot(fusion_analyzer_data(:, 25), '-m'); hold on;
    h(3) = plot(fusion_analyzer_data(:, 26), '-y');
    set(h,'LineWidth',4);
    grid on; legend('Lrate w.r.t S2','Lrate w.r.t. R1', 'Lrate w.r.t. R2');
    % ------------------m3-----------------
    heta3 = subplot(2,2,3);
    h(1) = plot(fusion_analyzer_data(:, 27), '-k'); hold on;
    h(2) = plot(fusion_analyzer_data(:, 28), '-y');
    set(h,'LineWidth',4);
    grid on; legend('Lrate w.r.t S3','Lrate w.r.t. R2');
    % ------------------m4-----------------
    heta4 = subplot(2,2,4);
    h(1) = plot(fusion_analyzer_data(:, 29), '-k'); hold on;
    h(2) = plot(fusion_analyzer_data(:, 30), '-y');
    set(h,'LineWidth',4);
    grid on; legend('Lrate w.r.t S4','Lrate w.r.t. R2' );
    
    heta = [heta1 heta2 heta3 heta4];
    % link axes for analysis
    linkaxes(heta, 'x');
    
    % set figure props
    set(gcf,'color','w');
end

%% plot the error gradients
if(strcmp(confidence_factor_type,'complex')==1)
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
figure;
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

%% COMPARATIVE VISUALIZATION FOR CONFIDENCE FACTOR ADAPTATION ANALYSIS
% load datasets
fixed = load('network_dataset_fixed.mat');
decay = load('network_dataset_decay.mat');
adaptive = load('network_dataset_adaptive.mat');
divisive = load('network_dataset_divisive.mat');
grad_hist = load('network_dataset_grad-history.mat');
grad_hist_avg = load('network_dataset_grad-history-average.mat');

% ---------m1----------
figure;
subplot(1,2,1);
plot(fixed.net_data(:,7), '.r'); hold on;
plot(decay.net_data(:,7),'.g'); hold on;
plot(divisive.net_data(:,7),'.b'); hold on;
plot(adaptive.net_data(:,7),'.k'); hold on;
plot(grad_hist.net_data(:,7),'.m'); hold on;
plot(grad_hist_avg.net_data(:,7),'.y');
title('Err w.r.t. S');
legend('Fixed cf', 'Decay cf', 'Divisive cf', 'Adaptive cf', 'Grad Hist cf', 'Grad Hist Avg cf');
grid on;
subplot(1,2,2);
plot(fixed.net_data(:,8), '.r'); hold on;
plot(decay.net_data(:,8),'.g'); hold on;
plot(divisive.net_data(:,8),'.b'); hold on;
plot(adaptive.net_data(:,8),'.k'); hold on;
plot(grad_hist.net_data(:,8),'.m'); hold on;
plot(grad_hist_avg.net_data(:,8),'.y');
title('Err w.r.t. R1');
legend('Fixed cf', 'Decay cf', 'Divisive cf', 'Adaptive cf', 'Grad Hist cf', 'Grad Hist Avg cf');
grid on;
% set figure props
set(gcf,'color','w');

% ---------m2----------
figure;
subplot(1,3,1);
plot(fixed.net_data(:,9), '.r'); hold on;
plot(decay.net_data(:,9),'.g'); hold on;
plot(divisive.net_data(:,9),'.b'); hold on;
plot(adaptive.net_data(:,9),'.k'); hold on;
plot(grad_hist.net_data(:,9),'.m'); hold on;
plot(grad_hist_avg.net_data(:,9),'.y');
title('Err w.r.t. S');
legend('Fixed cf', 'Decay cf', 'Divisive cf', 'Adaptive cf', 'Grad Hist cf', 'Grad Hist Avg cf');
grid on;

subplot(1,3,2);
plot(fixed.net_data(:,10), '.r'); hold on;
plot(decay.net_data(:,10),'.g'); hold on;
plot(divisive.net_data(:,10),'.b'); hold on;
plot(adaptive.net_data(:,10),'.k'); hold on;
plot(grad_hist.net_data(:,10),'.m'); hold on;
plot(grad_hist_avg.net_data(:,10),'.y');
title('Err w.r.t. R1');
legend('Fixed cf', 'Decay cf', 'Divisive cf', 'Adaptive cf', 'Grad Hist cf', 'Grad Hist Avg cf');
grid on;

subplot(1,3,3);
plot(fixed.net_data(:,11), '.r'); hold on;
plot(decay.net_data(:,11),'.g'); hold on;
plot(divisive.net_data(:,11),'.b'); hold on;
plot(adaptive.net_data(:,11),'.k'); hold on;
plot(grad_hist.net_data(:,11),'.m'); hold on;
plot(grad_hist_avg.net_data(:,11),'.y');
title('Err w.r.t. R2');
legend('Fixed cf', 'Decay cf', 'Divisive cf', 'Adaptive cf', 'Grad Hist cf', 'Grad Hist Avg cf');
grid on;

% set figure props
set(gcf,'color','w');

% ---------m3----------
figure;
subplot(1,2,1);
plot(fixed.net_data(:,12), '.r'); hold on;
plot(decay.net_data(:,12),'.g'); hold on;
plot(divisive.net_data(:,12),'.b'); hold on;
plot(adaptive.net_data(:,12),'.k'); hold on;
plot(grad_hist.net_data(:,12),'.m'); hold on;
plot(grad_hist_avg.net_data(:,12),'.y');
title('Err w.r.t. S');
legend('Fixed cf', 'Decay cf', 'Divisive cf', 'Adaptive cf', 'Grad Hist cf', 'Grad Hist Avg cf');
grid on;
subplot(1,2,2);
plot(fixed.net_data(:,13), '.r'); hold on;
plot(decay.net_data(:,13),'.g'); hold on;
plot(divisive.net_data(:,13),'.b'); hold on;
plot(adaptive.net_data(:,13),'.k'); hold on;
plot(grad_hist.net_data(:,13),'.m'); hold on;
plot(grad_hist_avg.net_data(:,13),'.y');
title('Err w.r.t. R2');
legend('Fixed cf', 'Decay cf', 'Divisive cf', 'Adaptive cf', 'Grad Hist cf', 'Grad Hist Avg cf');
grid on;
% set figure props
set(gcf,'color','w');

% ---------m4----------
figure;
subplot(1,2,1);
plot(fixed.net_data(:,14), '.r'); hold on;
plot(decay.net_data(:,14),'.g'); hold on;
plot(divisive.net_data(:,14),'.b'); hold on;
plot(adaptive.net_data(:,14),'.k'); hold on;
plot(grad_hist.net_data(:,14),'.m'); hold on;
plot(grad_hist_avg.net_data(:,14),'.y');
title('Err w.r.t. S');
legend('Fixed cf', 'Decay cf', 'Divisive cf', 'Adaptive cf', 'Grad Hist cf', 'Grad Hist Avg cf');
grid on;
subplot(1,2,2);
plot(fixed.net_data(:,15), '.r'); hold on;
plot(decay.net_data(:,15),'.g'); hold on;
plot(divisive.net_data(:,15),'.b'); hold on;
plot(adaptive.net_data(:,15),'.k'); hold on;
plot(grad_hist.net_data(:,15),'.m'); hold on;
plot(grad_hist_avg.net_data(:,15),'.y');
title('Err w.r.t. R2');
legend('Fixed cf', 'Decay cf', 'Divisive cf', 'Adaptive cf', 'Grad Hist cf', 'Grad Hist Avg cf');
grid on;
% set figure props
set(gcf,'color','w');