%% Fusion network for testing simple algebraic and temporal relationships
close all;
clear all;
clc;

%% Initialization

% enable or disable visualization
visualization_on = 1;

% prepare sensor data  to fed to the net
run_steps = 5000;
show_step  = 1000;
run_iters = run_steps - 1;
run_iters_extended = 2*run_iters;
maps = 6;                           % number of maps

% sensor values init - each sensor has it's own model
m1_sensor = build_data_set(0.523, show_step, run_iters);
m2_sensor = build_data_set(0.404, show_step, run_iters);
m3_sensor = build_data_set(0.627, show_step, run_iters);
m4_sensor = build_data_set(0.834, show_step, run_iters);
m5_sensor = build_data_set(0.702, show_step, run_iters);
m6_sensor = build_data_set(0.954, show_step, run_iters);

% complex / simple learning rate adaptation rules
complex_learning_update = 'simple'; % {simple, complex}
% params
switch complex_learning_update
    case 'simple'
        err_type = 'none';
        type = 'decay';                 % {decay, divisive}
    case 'complex'
        err_type = 'squared';            % {squared}
        type = 'up-down-factor';        % {up-down-factor, min-max-factor}
end

rand_update_rule = 1;           % id of the current update rule
update_rules = 14;              % possible update rules - depends on net topology
rules_ids = 1:update_rules;     % individual id for each rule
convergence_steps = 1;          % convergence steps counter
sim_points = run_iters_extended;

% sensor connection state - connected to net
sensor_connected = ones(1, maps);

% network iterator
net_iter = 1;

% init maps and indices
m1 = rand; m2 = rand; m3 = rand;
m4 = rand; m5 = rand; m6 = rand;

% maps ids
m1_id = 1; m2_id = 2; m3_id = 3;
m4_id = 4; m5_id = 5; m6_id = 6;

% number of connections of the maps (sensor + relations to other maps)
m1_links = 1:2; m2_links = 1:3;
m3_links = 1:2; m4_links = 1:3;
m5_links = 1:2; m6_links = 1:2;

% init maps errors for learning rate adaptation
em1 = zeros(1, length(m1_links)); em2 = zeros(1, length(m2_links));
em3 = zeros(1, length(m3_links)); em4 = zeros(1, length(m4_links));
em5 = zeros(1, length(m5_links)); em6 = zeros(1, length(m6_links));

% error gradients w.r.t. each map initialization (for learning rate adaptation)
dem1 = zeros(1, length(m1_links)); dem2 = zeros(1, length(m2_links));
dem3 = zeros(1, length(m3_links)); dem4 = zeros(1, length(m4_links));
dem5 = zeros(1, length(m5_links)); dem6 = zeros(1, length(m6_links));

dem1_old = zeros(1, length(m1_links)); dem2_old = zeros(1, length(m2_links));
dem3_old = zeros(1, length(m3_links)); dem4_old = zeros(1, length(m4_links));
dem5_old = zeros(1, length(m5_links)); dem6_old = zeros(1, length(m6_links));

% net params
% learning rates setup
maps_nr = maps;
sensors_nr = maps_nr;
error_nr = update_rules;
lrates_nr = update_rules;
net_data = zeros(sim_points, maps_nr+error_nr+lrates_nr);

% init learning rates
% bounds
ETA = 0.002;
ETAH = 10*ETA;
% init
etam1 = ETA*ones(sim_points, length(m1_links));
etam2 = ETA*ones(sim_points, length(m2_links));
etam3 = ETA*ones(sim_points, length(m3_links));
etam4 = ETA*ones(sim_points, length(m4_links));
etam5 = ETA*ones(sim_points, length(m5_links));
etam6 = ETA*ones(sim_points, length(m6_links));

%% Network dynamics
% encoded relationships
%
%   m2 = 3*m1
%   m3 = m2*m4
%   m4 = m5 + 2*m6

while(1)
    % decouple the sensors to check net relaxation
    if (convergence_steps == run_iters+1)
        sensor_connected = zeros(1, maps);
    end;
    % run the net to prove convergence after sensor decoupling
    if (convergence_steps ==  run_iters_extended + 1)
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
    % go through current sequence
    for i = 1:update_rules
        % pick a rule to update
        rand_update_rule = rules_ids(i);
        
        % update the corresponding node
        switch rand_update_rule
            %% m1 updates
            case 1
                if(sensor_connected(m1_id)==1)
                    m1 = (1-ETA)*m1 + etam1(convergence_steps, m1_links(1))*...
                        m1_sensor(net_iter);
                end
            case 2
                m1 = (1-ETA)*m1 + etam1(convergence_steps, m1_links(2)) * ...
                    1/3*m2;
                %% m2 updates
            case 3
                if(sensor_connected(m2_id)==1)
                    m2 = (1-ETA)*m2 + etam2(convergence_steps, m2_links(1))*...
                        m2_sensor(net_iter);
                end
            case 4
                m2 = (1-ETA)*m2 + etam2(convergence_steps, m2_links(2)) * ...
                    3*m1;
            case 5
                m2 = (1-ETA)*m2 + etam2(convergence_steps, m2_links(3)) * ...
                    m3/m4;
                %% m3 updates
            case 6
                if(sensor_connected(m3_id)==1)
                    m3 = (1-ETA)*m3 + etam3(convergence_steps, m3_links(1))*...
                        m3_sensor(net_iter);
                end
            case 7
                m3 = (1-ETA)*m3 + etam3(convergence_steps, m3_links(2)) * ...
                    m2*m4;
                %% m4 updates
            case 8
                if(sensor_connected(m4_id)==1)
                    m4 = (1-ETA)*m4 + etam4(convergence_steps, m4_links(1))*...
                        m4_sensor(net_iter);
                end
            case 9
                m4 = (1-ETA)*m4 + etam4(convergence_steps, m4_links(2)) * ...
                    m3/m2;
            case 10
                m4 = (1-ETA)*m4 + etam4(convergence_steps, m4_links(3)) * ...
                    (m5 + 2*m6);
                %% m5 updates
            case 11
                if(sensor_connected(m5_id)==1)
                    m5 = (1-ETA)*m5 + etam5(convergence_steps, m5_links(1))*...
                        m5_sensor(net_iter);
                end
            case 12
                m5 = (1-ETA)*m5 + etam5(convergence_steps, m5_links(2)) * ...
                    (m4 - 2*m6);
                %% m6 updates
            case 13
                if(sensor_connected(m6_id)==1)
                    m6 = (1-ETA)*m6 + etam6(convergence_steps, m6_links(1))*...
                        m6_sensor(net_iter);
                end
            case 14
                m6 = (1-ETA)*m6 + etam6(convergence_steps, m6_links(2)) * ...
                    ((m4 - m5)/2);
        end                
    end
    
                    % --------errors--------
                % error computation for learning rate adaptation
                % m1
                if(sensor_connected(m1_id)==1)
                    em1(1) = m1 - m1_sensor(net_iter);
                end
                em1(2) = m1 - m2/3;
                % m2
                if(sensor_connected(m2_id)==1)
                    em2(1) = m2 - m2_sensor(net_iter);
                end
                em2(2) = m2 - 3*m1;
                em2(3) = m2 - m3/m4;
                % m3
                if(sensor_connected(m3_id)==1)
                    em3(1) = m3 - m3_sensor(net_iter);
                end
                em3(2) = m3 - m2*m4;
                % m4
                if(sensor_connected(m4_id)==1)
                    em4(1) = m4 - m4_sensor(net_iter);
                end
                em4(2) = m4 - m3/m2;
                em4(3) = m4 - (m5+2*m6);
                % m5
                if(sensor_connected(m5_id)==1)
                    em5(1) = m5 - m5_sensor(net_iter);
                end
                em5(2) = m5 - (m4 - 2*m6);
                % m6
                if(sensor_connected(m6_id)==1)
                    em6(1) = m6 - m6_sensor(net_iter);
                end
                em6(2) = m6 - ((m4 - m5)/2);
                
            %% learning rates adaptation and clamping
                    % compute the error gradients for simple / squared error for
        % complex update rules
        switch err_type
              
            case 'squared'
                
                % -----------errors------------
                % error computation for learning rate adaptation
                % m1
                if(sensor_connected(m1_id)==1)
                    em1(1) = (m1 - m1_sensor(net_iter))^2;
                end
                em1(2) = (m1 - m2/3)^2;
                % m2
                if(sensor_connected(m2_id)==1)
                    em2(1) = (m2 - m2_sensor(net_iter))^2;
                end
                em2(2) = (m2 - 3*m1)^2;
                em2(3) = (m2 - m3/m4)^2;
                % m3
                if(sensor_connected(m3_id)==1)
                    em3(1) = (m3 - m3_sensor(net_iter))^2;
                end
                em3(2) = (m3 - m2*m4)^2;
                % m4
                if(sensor_connected(m4_id)==1)
                    em4(1) = (m4 - m4_sensor(net_iter))^2;
                end
                em4(2) = (m4 - m3/m2)^2;
                em4(3) = (m4 - (m5+2*m6))^2;
                % m5
                if(sensor_connected(m5_id)==1)
                    em5(1) = (m5 - m5_sensor(net_iter))^2;
                end
                em5(2) = (m5 - (m4 - 2*m6))^2;
                % m6
                if(sensor_connected(m6_id)==1)
                    em6(1) = (m6 - m6_sensor(net_iter))^2;
                end
                em6(2) = (m6 - ((m4 - m5)/2))^2;
                
                % ---------gradients----------
                % gradient of errors for m1
                if(sensor_connected(m1_id)==1)
                    dem1(m1_links(1)) = 2*(m1-m1_sensor(net_iter))*(-m1_sensor(net_iter));
                end
                dem1(m1_links(2)) = 2*(m1-m2/3)*(-m2/3);
                
                % gradient of errors for m2
                if(sensor_connected(m2_id)==1)
                    dem2(m2_links(1)) = 2*(m2-m2_sensor(net_iter))*(-m2_sensor(net_iter));
                end
                dem2(m2_links(2)) = 2*(m2-3*m1)*(-3*m1);
                dem2(m2_links(3)) = 2*(m2-m3/m4)*(-m3/m4);
                
                % gradient of errors for m3
                if(sensor_connected(m3_id)==1)
                    dem3(m3_links(1)) = 2*(m3-m3_sensor(net_iter))*(-m3_sensor(net_iter));
                end
                dem3(m3_links(2)) = 2*(m3-m2*m4)*(-m2/m4);
                
                % gradient of errors for m4
                if(sensor_connected(m4_id)==1)
                    dem4(m4_links(1)) = 2*(m4-m4_sensor(net_iter))*(-m4_sensor(net_iter));
                end
                dem4(m4_links(2)) = 2*(m4-m3/m2)*(-m3/m2);
                dem4(m4_links(3)) = 2*(m4-(m5+2*m6))*(-(m5+2*m6));
                
                % gradient of errors for m5
                if(sensor_connected(m5_id)==1)
                    dem5(m5_links(1)) = 2*(m5-m5_sensor(net_iter))*(-m5_sensor(net_iter));
                end
                dem5(m5_links(2)) = 2*(m5-(m4-2*m6))*(-(m4-2*m6));
                
                % gradient of errors for m6
                if(sensor_connected(m6_id)==1)
                    dem6(m6_links(1)) = 2*(m6-m6_sensor(net_iter))*(-m6_sensor(net_iter));
                end
                dem6(m6_links(2)) = 2*(m6-((m4-m5)/2))*(-((m4-m5)/2));
        end
        
        switch(complex_learning_update)
            case 'simple'
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
                for k=1:length(m5_links)
                    etam5(convergence_steps, m5_links(k)) = update_learning_rate(etam5(convergence_steps, m5_links(k)), em5, m5_links(k), ETA, type);
                    etam5(convergence_steps, m5_links(k)) = clamp(etam5(convergence_steps, m5_links(k)), ETAH);
                end
                for k=1:length(m6_links)
                    etam6(convergence_steps, m6_links(k)) = update_learning_rate(etam6(convergence_steps, m6_links(k)), em6, m6_links(k), ETA, type);
                    etam6(convergence_steps, m6_links(k)) = clamp(etam6(convergence_steps, m6_links(k)), ETAH);
                end
            case 'complex'
                % params init
                u = 134;           % u > 1
                d = 0.00030554366;      % d < 1
                l_min = ETA;
                l_max = ETAH;
                
                % complex update rules for the learning rate
                for k=1:length(m1_links)
                    etam1(convergence_steps, m1_links(k)) = update_learning_rate_complex(etam1(convergence_steps, m1_links(k)), dem1(k), dem1_old(k), u, d, l_min, l_max, type);
                    etam1(convergence_steps, m1_links(k)) = clamp(etam1(convergence_steps, m1_links(k)), ETAH);
                end
                for k=1:length(m2_links)
                    etam2(convergence_steps, m2_links(k)) = update_learning_rate_complex(etam2(convergence_steps, m2_links(k)), dem2(k), dem2_old(k), u, d, l_min, l_max, type);
                    etam2(convergence_steps, m2_links(k)) = clamp(etam2(convergence_steps, m2_links(k)), ETAH);
                end
                
                for k=1:length(m3_links)
                    etam3(convergence_steps, m3_links(k)) = update_learning_rate_complex(etam3(convergence_steps, m3_links(k)), dem3(k), dem3_old(k), u, d, l_min, l_max, type);
                    etam3(convergence_steps, m3_links(k)) = clamp(etam3(convergence_steps, m3_links(k)), ETAH);
                end
                
                for k=1:length(m4_links)
                    etam4(convergence_steps, m4_links(k)) = update_learning_rate_complex(etam4(convergence_steps, m4_links(k)), dem4(k), dem4_old(k), u, d, l_min, l_max, type);
                    etam4(convergence_steps, m4_links(k)) = clamp(etam4(convergence_steps, m4_links(k)), ETAH);
                end
                
                for k=1:length(m5_links)
                    etam5(convergence_steps, m5_links(k)) = update_learning_rate_complex(etam5(convergence_steps, m5_links(k)), dem5(k), dem5_old(k), u, d, l_min, l_max, type);
                    etam5(convergence_steps, m5_links(k)) = clamp(etam5(convergence_steps, m5_links(k)), ETAH);
                end
                
                for k=1:length(m6_links)
                    etam6(convergence_steps, m6_links(k)) = update_learning_rate_complex(etam6(convergence_steps, m6_links(k)), dem6(k), dem6_old(k), u, d, l_min, l_max, type);
                    etam6(convergence_steps, m6_links(k)) = clamp(etam6(convergence_steps, m6_links(k)), ETAH);
                end
        end
        
        % update indices
        dem1_old = dem1; dem2_old = dem2; dm3_old = dem3;
        dem4_old = dem4; dem5_old = dem5; dm6_old = dem6;
        
    %% write data to net struct
    % maps
    net_data(convergence_steps, 1) = m1;
    net_data(convergence_steps, 2) = m2;
    net_data(convergence_steps, 3) = m3;
    net_data(convergence_steps, 4) = m4;
    net_data(convergence_steps, 5) = m5;
    net_data(convergence_steps, 6) = m6;
    % error
    net_data(convergence_steps, 7) = em1(1);
    net_data(convergence_steps, 8) = em1(2);
    
    net_data(convergence_steps, 9) = em2(1);
    net_data(convergence_steps, 10) = em2(2);
    net_data(convergence_steps, 11) = em2(3);
    
    net_data(convergence_steps, 12) = em3(1);
    net_data(convergence_steps, 13) = em3(2);
    
    net_data(convergence_steps, 14) = em4(1);
    net_data(convergence_steps, 15) = em4(2);
    net_data(convergence_steps, 16) = em4(3);
    
    net_data(convergence_steps, 17) = em5(1);
    net_data(convergence_steps, 18) = em5(2);
    
    net_data(convergence_steps, 19) = em6(1);
    net_data(convergence_steps, 20) = em6(2);
    
    % sample idx
    net_data(convergence_steps, 21) = convergence_steps;
    
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
    net_data(convergence_steps, 31) = etam4(convergence_steps, m4_links(3));
    
    net_data(convergence_steps, 32) = etam5(convergence_steps, m5_links(1));
    net_data(convergence_steps, 33) = etam5(convergence_steps, m5_links(2));
    
    net_data(convergence_steps, 34) = etam6(convergence_steps, m6_links(1));
    net_data(convergence_steps, 35) = etam6(convergence_steps, m6_links(2));
    
    %% update indices
    net_iter = net_iter + 1;
    convergence_steps = convergence_steps + 1;
end

%% fill in the net simulation data into the visualization struct
fusion_analyzer_data = net_data;

%% Visualization
if(visualization_on==1)
    figure(1);
    % ---------------- R1 ----------------------
    hd(1) = subplot(4, 3, 1);
    plot(net_data(:, 1), '.b');
    grid on; title('M1 map values');
    
    hd(2) = subplot(4, 3, 4);
    plot(net_data(:, 2), '.b');
    grid on; title('M2 map values');
    
    hd(3) = subplot(4, 3, 7);
    plot(m1_sensor, '.k');
    grid on; title('S1 sensor values');
    
    hd(4) = subplot(4, 3, 10);
    plot(m2_sensor, '.k');
    grid on; title('S2 sensor values');
    % -------------------------------------------
    
    % ---------------- R2 ----------------------
    hd(5) = subplot(4, 3, 2);
    plot(net_data(:, 3), '.r');
    grid on; title('M3 map values');
    
    hd(6) = subplot(4, 3, 5);
    plot(net_data(:, 4), '.r');
    grid on; title('M4 map values');
    
    hd(7) = subplot(4, 3, 8);
    plot(m3_sensor, '.k');
    grid on; title('S3 sensor values');
    
    hd(8) = subplot(4, 3, 11);
    plot(m2_sensor, '.k');
    grid on; title('S4 sensor values');
    % -------------------------------------------
    
    % ---------------- R3 ----------------------
    hd(9) = subplot(4, 3, 3);
    plot(net_data(:, 5), '.g');
    grid on; title('M5 map values');
    
    hd(10) = subplot(4, 3, 6);
    plot(net_data(:, 6), '.g');
    grid on; title('M6 map values');
    
    hd(11) = subplot(4, 3, 9);
    plot(m5_sensor, '.k');
    grid on; title('S5 sensor values');
    
    hd(12) = subplot(4, 3, 12);
    plot(m6_sensor, '.k');
    grid on; title('S6 sensor values');
    
    % link axes
    linkaxes(hd, 'x');
    
    % -------------- Erros signals --------------
    figure(2);
    %-------------------M1-------------------
    he(1) = subplot(2,3,1);
    plot(net_data(:,7),'.k');hold on;
    plot(net_data(:,8), '.b');
    title('M1 errors');
    legend('Err w.r.t S1', 'Err w.r.t R1');
    grid on;
    %-------------------M2-------------------
    he(2) = subplot(2,3,2);
    plot(net_data(:,9), '.k'); hold on;
    plot(net_data(:,10), '.b'); hold on;
    plot(net_data(:,11), '.r');
    title('M2 errors');
    legend('Err w.r.t S2', 'Err w.r.t R1', 'Err w.r.t R2');
    grid on;
    %-------------------M3-------------------
    he(3) = subplot(2,3,3);
    plot(net_data(:,12), '.k'); hold on;
    plot(net_data(:, 13), '.r');
    title('M3 errors');
    legend('Err w.r.t S3', 'Err w.r.t R2');
    grid on;
    %-------------------M4-------------------
    he(4) = subplot(2,3,4);
    plot(net_data(:,14), '.k'); hold on;
    plot(net_data(:,15), '.r'); hold on;
    plot(net_data(:, 16), '.g');
    title('M4 errors');
    legend('Err w.r.t S4', 'Err w.r.t R2', 'Err w.r.t R3');
    grid on;
    %-------------------M5-------------------
    he(5) = subplot(2,3,5);
    plot(net_data(:,17), '.k'); hold on;
    plot(net_data(:, 18), '.g');
    title('M5 errors');
    legend('Err w.r.t S5', 'Err w.r.t R3');
    grid on;
    %-------------------M6-------------------
    he(6) = subplot(2,3,6);
    plot(net_data(:,19), '.k'); hold on;
    plot(net_data(:, 20), '.g');
    title('M6 errors');
    legend('Err w.r.t S6', 'Err w.r.t R3');
    grid on;
    
    linkaxes(he, 'x');
    % -------------------------------------------
    
    %--------------------------------------------------------------------------------------------------------
    fig2Handle = figure(3);
    set(fig2Handle, 'Position', [100, 100, 600, 1000]);
    %-------------------------
    subplot(3,1,1)
    % axis control and adjustment
    [ax1] = axescheck(fusion_analyzer_data(:,20), fusion_analyzer_data(:,1));
    ax1 = newplot(ax1);
    set(ax1,'XGrid','on');
    set(ax1,'YGrid','on');
    if ~ishold(ax1)
        [minx1,maxx1] = minmax(fusion_analyzer_data(:,20));
        [miny1,maxy1] = minmax(fusion_analyzer_data(:,1));
        axis(ax1,[minx1 maxx1 miny1 maxy1])
    end
    title('Maps values during simulation');
    ylabel('Map 1 values')
    xlabel('Data points')
    hc1 =  line('parent',ax1, 'linestyle',':','LineWidth', 2.0,'erase','xor','xdata',fusion_analyzer_data(1,20),'ydata',fusion_analyzer_data(1,1));
    ha1 = line('parent',ax1,'linestyle',':','LineWidth', 2.0,'erase','none','xdata',[],'ydata',[]);
    %-------------------------
    subplot(3,1,2)
    % axis control and adjustment
    [ax2] = axescheck(fusion_analyzer_data(:,20), fusion_analyzer_data(:,2));
    ax2 = newplot(ax2);
    set(ax2,'XGrid','on');
    set(ax2,'YGrid','on');
    if ~ishold(ax2)
        [minx2,maxx2] = minmax(fusion_analyzer_data(:,20));
        [miny2,maxy2] = minmax(fusion_analyzer_data(:,2));
        axis(ax2,[minx2 maxx2 miny2 maxy2])
    end
    ylabel('Map 2 values')
    xlabel('Data points')
    hc2 =  line('parent',ax2,'linestyle',':','LineWidth', 2.0,'erase','xor','xdata',fusion_analyzer_data(1,20),'ydata',fusion_analyzer_data(1,2));
    ha2 = line('parent',ax2,'linestyle',':','LineWidth', 2.0,'erase','none','xdata',[],'ydata',[]);
    %-------------------------
    subplot(3,1,3)
    % axis control and adjustment
    [ax3] = axescheck(fusion_analyzer_data(:,1), fusion_analyzer_data(:,2));
    ax3 = newplot(ax3);
    if ~ishold(ax3)
        [minx3,maxx3] = minmax(fusion_analyzer_data(:,1));
        [miny3,maxy3] = minmax(fusion_analyzer_data(:,2));
        axis(ax3,[minx3 maxx3 miny3 maxy3])
    end
    % interval for the expected trajectory
    abs_max = max(abs(minx3), abs(maxx3));
    t = -abs_max:abs_max;
    plot(t, 3*t, '-m');
    set(ax3,'XGrid','on');
    set(ax3,'YGrid','on');
    hc3 =  line('parent',ax3,'linestyle',':','LineWidth', 2.0,'erase','xor','xdata',fusion_analyzer_data(1,1),'ydata',fusion_analyzer_data(1,2));
    ha3 = line('parent',ax3,'linestyle',':','LineWidth', 2.0,'erase','none','xdata',[],'ydata',[]);
    title('M2 dependency on M1');
    xlabel('M1 values');
    ylabel('M2 values');
    % dynamic plot simulataneously in all subplots
    for i=2:length(fusion_analyzer_data(:,20))
        j = i-1:i;
        set(hc1,'xdata', fusion_analyzer_data(i,20), 'ydata', fusion_analyzer_data(i,1));
        set(ha1,'xdata', fusion_analyzer_data(j,20), 'ydata', fusion_analyzer_data(j,1));
        set(hc2,'xdata', fusion_analyzer_data(i,20), 'ydata', fusion_analyzer_data(i,2));
        set(ha2,'xdata', fusion_analyzer_data(j,20), 'ydata', fusion_analyzer_data(j,2));
        set(hc3,'xdata', fusion_analyzer_data(i,1), 'ydata', fusion_analyzer_data(i,2));
        set(ha3,'xdata', fusion_analyzer_data(j,1), 'ydata', fusion_analyzer_data(j,2));
        drawnow;
    end
    
    %--------------------------------------------------------------------------------------------------------
    fig3Handle = figure(4);
    set(fig3Handle, 'Position', [900, 900, 600, 1200]);
    subplot(4,1,1)
    % axis control and adjustment
    [ax4] = axescheck(fusion_analyzer_data(:,20), fusion_analyzer_data(:,3));
    ax4 = newplot(ax4);
    set(ax4,'XGrid','on');
    set(ax4,'YGrid','on');
    if ~ishold(ax4)
        [minx4,maxx4] = minmax(fusion_analyzer_data(:,20));
        [miny4,maxy4] = minmax(fusion_analyzer_data(:,3));
        axis(ax4,[minx4 maxx4 miny4 maxy4])
    end
    title('Maps values during simulation');
    hc4 =  line('parent',ax4, 'linestyle',':','LineWidth', 2.0,'erase','xor','xdata',fusion_analyzer_data(1,20),'ydata',fusion_analyzer_data(1,3));
    ha4 = line('parent',ax4,'linestyle',':','LineWidth', 2.0,'erase','none','xdata',[],'ydata',[]);
    ylabel('Map 3 values')
    xlabel('Data points')
    
    %-------------
    subplot(4,1,2)
    % axis control and adjustment
    [ax5] = axescheck(fusion_analyzer_data(:,20), fusion_analyzer_data(:,4));
    ax5 = newplot(ax5);
    set(ax5,'XGrid','on');
    set(ax5,'YGrid','on');
    if ~ishold(ax5)
        [minx5,maxx5] = minmax(fusion_analyzer_data(:,20));
        [miny5,maxy5] = minmax(fusion_analyzer_data(:,4));
        axis(ax5,[minx5 maxx5 miny5 maxy5])
    end
    hc5 =  line('parent',ax5, 'linestyle',':','LineWidth', 2.0,'erase','xor','xdata',fusion_analyzer_data(1,20),'ydata',fusion_analyzer_data(1,4));
    ha5 = line('parent',ax5,'linestyle',':','LineWidth', 2.0,'erase','none','xdata',[],'ydata',[]);
    ylabel('Map 4 values')
    xlabel('Data points')
    %-------------
    subplot(4,1,3)
    % axis control and adjustment
    [ax6] = axescheck(fusion_analyzer_data(:,20), fusion_analyzer_data(:,5));
    ax6 = newplot(ax6);
    set(ax6,'XGrid','on');
    set(ax6,'YGrid','on');
    if ~ishold(ax6)
        [minx6,maxx6] = minmax(fusion_analyzer_data(:,20));
        [miny6,maxy6] = minmax(fusion_analyzer_data(:,5));
        axis(ax6,[minx6 maxx6 miny6 maxy6])
    end
    hc6 = line('parent',ax6, 'linestyle',':','LineWidth', 2.0,'erase','xor','xdata',fusion_analyzer_data(1,20),'ydata',fusion_analyzer_data(1,5));
    ha6 = line('parent',ax6,'linestyle',':','LineWidth', 2.0,'erase','none','xdata',[],'ydata',[]);
    ylabel('Map 5 values')
    xlabel('Data points')
    %------------
    subplot(4,1,4)
    % axis control and adjustment
    [ax7] = axescheck(fusion_analyzer_data(:,20), fusion_analyzer_data(:,6));
    ax7 = newplot(ax7);
    set(ax7,'XGrid','on');
    set(ax7,'YGrid','on');
    if ~ishold(ax7)
        [minx7,maxx7] = minmax(fusion_analyzer_data(:,20));
        [miny7,maxy7] = minmax(fusion_analyzer_data(:,6));
        axis(ax7,[minx7 maxx7 miny7 maxy7])
    end
    hc7 =  line('parent',ax7, 'linestyle',':','LineWidth', 2.0,'erase','xor','xdata',fusion_analyzer_data(1,20),'ydata',fusion_analyzer_data(1,6));
    ha7 = line('parent',ax7,'linestyle',':','LineWidth', 2.0,'erase','none','xdata',[],'ydata',[]);
    ylabel('Map 6 values')
    xlabel('Data points')
    % loop and dynamically shouw map convergence
    for i=2:length(fusion_analyzer_data(:,20))
        j = i-1:i;
        set(hc4,'xdata', fusion_analyzer_data(i,20), 'ydata', fusion_analyzer_data(i,3));
        set(ha4,'xdata', fusion_analyzer_data(j,20), 'ydata', fusion_analyzer_data(j,3));
        set(hc5,'xdata', fusion_analyzer_data(i,20), 'ydata', fusion_analyzer_data(i,4));
        set(ha5,'xdata', fusion_analyzer_data(j,20), 'ydata', fusion_analyzer_data(j,4));
        set(hc6,'xdata', fusion_analyzer_data(i,20), 'ydata', fusion_analyzer_data(i,5));
        set(ha6,'xdata', fusion_analyzer_data(j,20), 'ydata', fusion_analyzer_data(j,5));
        set(hc7,'xdata', fusion_analyzer_data(i,20), 'ydata', fusion_analyzer_data(i,6));
        set(ha7,'xdata', fusion_analyzer_data(j,20), 'ydata', fusion_analyzer_data(j,6));
        drawnow;
    end
    
    % %--------------------------------------------------------------------------------------------------------
    % % per relation errors
    % % first relation between M1 and M2
    % figure(5);
    % subplot(2,1,1)
    % plot(fusion_analyzer_data(:,20), fusion_analyzer_data(:,8))
    % ylabel('Map 1 error to R1')
    % xlabel('Data points')
    % grid on
    % subplot(2,1,2);
    % plot(fusion_analyzer_data(:,20), fusion_analyzer_data(:,10));
    % ylabel('Map 2 error to R1');
    % xlabel('Data points');
    % grid on
    % %--------------------------------------------------------------------------------------------------------
    % % second relationship between M2, M3 and M4
    % figure(6);
    % subplot(3,1,1)
    % plot(fusion_analyzer_data(:,20), fusion_analyzer_data(:,11))
    % ylabel('Map 2 error to R2');
    % xlabel('Data points');
    % grid on;
    % subplot(3,1,2)
    % plot(fusion_analyzer_data(:,20), fusion_analyzer_data(:,13))
    % ylabel('Map 3 error to R2');
    % xlabel('Data points');
    % grid on;
    % subplot(3,1,3)
    % plot(fusion_analyzer_data(:,20), fusion_analyzer_data(:,15))
    % ylabel('Map 4 error to R2');
    % xlabel('Data points');
    % grid on;
    % %--------------------------------------------------------------------------------------------------------
    % % third relationship between M4, M5, M6
    % figure(7);
    % subplot(3,1,1)
    % plot(fusion_analyzer_data(:,20), fusion_analyzer_data(:,16))
    % ylabel('Map 4 error to R3')
    % xlabel('Data points')
    % grid on
    % subplot(3,1,2)
    % plot(fusion_analyzer_data(:,20), fusion_analyzer_data(:,18))
    % ylabel('Map 5 error to R3')
    % xlabel('Data points')
    % grid on
    % subplot(3,1,3)
    % plot(fusion_analyzer_data(:,20), fusion_analyzer_data(:,20))
    % ylabel('Map 6 error to R3')
    % xlabel('Data points')
    % grid on
    
    %--------------------------------------------------------------------------------------------------------
    % analize the first relationship R1: M2 = 3*M1
    % analize the (M1, M2) dependency
    figure(8);
    plot(fusion_analyzer_data(:,20), 3*fusion_analyzer_data(:,1),'--r');
    hold on;
    plot(fusion_analyzer_data(:,20), fusion_analyzer_data(:,2),'--b');
    title('The dependency between M1 and M2');
    grid on;
    legend('3*M1 data', 'M2 data');
    xlabel('Data points');
    
    [ax8] = axescheck(fusion_analyzer_data(:,20), 3*fusion_analyzer_data(:,1));
    ax8 = newplot(ax8);
    set(ax8,'XGrid','on');
    set(ax8,'YGrid','on');
    if ~ishold(ax8)
        [minx8,maxx8] = minmax(fusion_analyzer_data(:,20));
        [miny8,maxy8] = minmax(3*fusion_analyzer_data(:,1));
        axis(ax8,[minx8 maxx8 miny8 maxy8])
    end
    hc8 = line('parent',ax8, 'linestyle',':','LineWidth', 2.0,'erase','xor','xdata',fusion_analyzer_data(1,20),'ydata',3*fusion_analyzer_data(1,1));
    ha8 = line('parent',ax8,'linestyle',':','LineWidth', 2.0,'erase','none','xdata',[],'ydata',[]);
    
    %------------
    % axis control and adjustment
    [ax9] = axescheck(fusion_analyzer_data(:,20), fusion_analyzer_data(:,2));
    ax9 = newplot(ax9);
    set(ax9,'XGrid','on');
    set(ax9,'YGrid','on');
    if ~ishold(ax9)
        [minx9,maxx9] = minmax(fusion_analyzer_data(:,20));
        [miny9,maxy9] = minmax(fusion_analyzer_data(:,2));
        axis(ax9,[minx9 maxx9 miny9 maxy9])
    end
    hc9 =  line('parent',ax9, 'linestyle',':','LineWidth', 2.0,'erase','xor','xdata',fusion_analyzer_data(1,20),'ydata',fusion_analyzer_data(1,2));
    ha9 = line('parent',ax9,'linestyle',':','LineWidth', 2.0,'erase','none','xdata',[],'ydata',[]);
    
    % loop and dynamically shouw map convergence
    for i=2:length(fusion_analyzer_data(:,20))
        j = i-1:i;
        set(hc8,'xdata', fusion_analyzer_data(i,20), 'ydata', 3*fusion_analyzer_data(i,1));
        set(ha8,'xdata', fusion_analyzer_data(j,20), 'ydata', 3*fusion_analyzer_data(j,1));
        set(hc9,'xdata', fusion_analyzer_data(i,20), 'ydata', fusion_analyzer_data(i,2));
        set(ha9,'xdata', fusion_analyzer_data(j,20), 'ydata', fusion_analyzer_data(j,2));
        drawnow;
    end
    
    %--------------------------------------------------------------------------------------------------------
    % analize the second relationship R2: M3 = M2*M4
    % analize the (M2, M3, M4) dependency
    figure(9);
    plot(fusion_analyzer_data(:,20), fusion_analyzer_data(:,2).*fusion_analyzer_data(:,4), '--r');
    hold on;
    plot(fusion_analyzer_data(:,20), fusion_analyzer_data(:,3), '--b');
    title('The dependency between M2, M3, and M4');
    legend('M2*M4 data', 'M3 data');
    grid on;
    xlabel('Data points');
    
    [ax10] = axescheck(fusion_analyzer_data(:,20), fusion_analyzer_data(:,2).*fusion_analyzer_data(:,4));
    ax10 = newplot(ax10);
    set(ax10,'XGrid','on');
    set(ax10,'YGrid','on');
    if ~ishold(ax10)
        [minx10,maxx10] = minmax(fusion_analyzer_data(:,20));
        [miny10,maxx10] = minmax(fusion_analyzer_data(:,2).*fusion_analyzer_data(:,4));
        axis(ax10,[minx10 maxx10 miny10 maxy10])
    end
    hc10 = line('parent',ax10, 'linestyle',':','LineWidth', 2.0,'erase','xor','xdata',fusion_analyzer_data(1,20),'ydata',fusion_analyzer_data(1,2).*fusion_analyzer_data(1,4));
    ha10 = line('parent',ax10,'linestyle',':','LineWidth', 2.0,'erase','none','xdata',[],'ydata',[]);
    
    %------------
    % axis control and adjustment
    [ax11] = axescheck(fusion_analyzer_data(:,20), fusion_analyzer_data(:,3));
    ax11 = newplot(ax11);
    set(ax11,'XGrid','on');
    set(ax11,'YGrid','on');
    if ~ishold(ax11)
        [minx11,maxx11] = minmax(fusion_analyzer_data(:,20));
        [miny11,maxy11] = minmax(fusion_analyzer_data(:,2));
        axis(ax11,[minx11 maxx11 miny11 maxy11])
    end
    hc11 =  line('parent',ax11, 'linestyle',':','LineWidth', 2.0,'erase','xor','xdata',fusion_analyzer_data(1,20),'ydata', fusion_analyzer_data(1,3));
    ha11 = line('parent',ax11,'linestyle',':','LineWidth', 2.0,'erase','none','xdata',[],'ydata',[]);
    
    % loop and dynamically shouw map convergence
    for i=2:length(fusion_analyzer_data(:,20))
        j = i-1:i;
        set(hc10,'xdata', fusion_analyzer_data(i,20), 'ydata', fusion_analyzer_data(i,2).*fusion_analyzer_data(i,4));
        set(ha10,'xdata', fusion_analyzer_data(j,20), 'ydata', fusion_analyzer_data(j,2).*fusion_analyzer_data(j,4));
        set(hc11,'xdata', fusion_analyzer_data(i,20), 'ydata', fusion_analyzer_data(i,3));
        set(ha11,'xdata', fusion_analyzer_data(j,20), 'ydata', fusion_analyzer_data(j,3));
        drawnow;
    end
    
    %--------------------------------------------------------------------------------------------------------
    % analize the second relationship R3: M4 = M5 + 2*M6
    % analize the (M4, M5, M6) dependency
    figure(10);
    plot(fusion_analyzer_data(:,20), fusion_analyzer_data(:,5)+(2*fusion_analyzer_data(:,6)), '--r');
    hold on;
    plot(fusion_analyzer_data(:,20), fusion_analyzer_data(:,4),'--b');
    title('The dependency between M4, M5, and M6');
    legend('M5+2*M6 data', 'M4 data');
    grid on;
    xlabel('Data points');
    
    [ax12] = axescheck(fusion_analyzer_data(:,20), fusion_analyzer_data(:,5)+ 2*fusion_analyzer_data(:,6));
    ax12 = newplot(ax12);
    set(ax12,'XGrid','on');
    set(ax12,'YGrid','on');
    if ~ishold(ax12)
        [minx12,maxx12] = minmax(fusion_analyzer_data(:,20));
        [miny12,maxx12] = minmax(fusion_analyzer_data(:,5)+2*fusion_analyzer_data(:,6));
        axis(ax12,[minx12 maxx12 miny12 maxy12])
    end
    hc12 = line('parent',ax12, 'linestyle',':','LineWidth', 2.0,'erase','xor','xdata',fusion_analyzer_data(1,20),'ydata',fusion_analyzer_data(1,5)+2*fusion_analyzer_data(1,6));
    ha12 = line('parent',ax12,'linestyle',':','LineWidth', 2.0,'erase','none','xdata',[],'ydata',[]);
    
    %------------
    % axis control and adjustment
    [ax13] = axescheck(fusion_analyzer_data(:,20), fusion_analyzer_data(:,4));
    ax13 = newplot(ax13);
    set(ax13,'XGrid','on');
    set(ax13,'YGrid','on');
    if ~ishold(ax13)
        [minx13,maxx13] = minmax(fusion_analyzer_data(:,20));
        [miny13,maxy13] = minmax(fusion_analyzer_data(:,4));
        axis(ax13,[minx13 maxx13 miny13 maxy13])
    end
    hc13 =  line('parent',ax13, 'linestyle',':','LineWidth', 2.0,'erase','xor','xdata',fusion_analyzer_data(1,20),'ydata', fusion_analyzer_data(1,4));
    ha13 = line('parent',ax13,'linestyle',':','LineWidth', 2.0,'erase','none','xdata',[],'ydata',[]);
    
    % loop and dynamically shouw map convergence
    for i=2:length(fusion_analyzer_data(:,20))
        j = i-1:i;
        set(hc12,'xdata', fusion_analyzer_data(i,20), 'ydata', fusion_analyzer_data(i,5)+2*fusion_analyzer_data(i,6));
        set(ha12,'xdata', fusion_analyzer_data(j,20), 'ydata', fusion_analyzer_data(j,5)+2*fusion_analyzer_data(j,6));
        set(hc13,'xdata', fusion_analyzer_data(i,20), 'ydata', fusion_analyzer_data(i,4));
        set(ha13,'xdata', fusion_analyzer_data(j,20), 'ydata', fusion_analyzer_data(j,4));
        drawnow;
    end
end % visualization on flag

% plot learning rates on a per map basis
figure(11);
% ------------------m1-----------------
heta1 = subplot(2,3,1);
plot(fusion_analyzer_data(:, 22), '.k'); hold on;
plot(fusion_analyzer_data(:, 23), '.r');
grid on;
legend('Lrate w.r.t S1','Lrate w.r.t. R1');
title('Learning rate analysis');
% ------------------m1-----------------
heta2 = subplot(2,3,2);
plot(fusion_analyzer_data(:, 24), '.k'); hold on;
plot(fusion_analyzer_data(:, 25), '.r'); hold on;
plot(fusion_analyzer_data(:, 26), '.b');
grid on; legend('Lrate w.r.t S2','Lrate w.r.t. R1', 'Lrate w.r.t. R2');
% ------------------m1-----------------
heta3 = subplot(2,3,3);
plot(fusion_analyzer_data(:, 27), '.k'); hold on;
plot(fusion_analyzer_data(:, 28), '.b');
grid on; legend('Lrate w.r.t S3','Lrate w.r.t. R2');
% ------------------m1-----------------
heta4 = subplot(2,3,4);
plot(fusion_analyzer_data(:, 29), '.k'); hold on;
plot(fusion_analyzer_data(:, 30), '.b'); hold on;
plot(fusion_analyzer_data(:, 31), '.g');
grid on; legend('Lrate w.r.t S4','Lrate w.r.t. R2', 'Lrate w.r.t. R3');
% ------------------m1-----------------
heta5 = subplot(2,3,5);
plot(fusion_analyzer_data(:, 32), '.k'); hold on;
plot(fusion_analyzer_data(:, 33), '.g');
grid on; legend('Lrate w.r.t S5','Lrate w.r.t. R3');
% ------------------m1-----------------
heta6 = subplot(2,3,6);
plot(fusion_analyzer_data(:, 34), '.k'); hold on;
plot(fusion_analyzer_data(:, 35), '.g');
grid on; legend('Lrate w.r.t S6', 'Lrate w.r.t. R3');

heta = [heta1 heta2 heta3 heta4 heta5 heta6];
% link axes for analysis
linkaxes(heta, 'x');

% %% save learning rate adaptation data in files for comparison
% f = fopen(strcat('eta_',type),'w');
% eta_all = [etam1 etam2 etam3 etam4 etam5 etam6];
% bytes = fprintf(f, '%6.6f %6.6f %6.6f %6.6f %6.6f %6.6f %6.6f %6.6f %6.6f %6.6f %6.6f %6.6f %6.6f %6.6f \n',eta_all);
% fclose(f);