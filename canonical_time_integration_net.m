%% Fusion network for testing simple integration and differentiation in time
% cleanup
close all;
clear all;
clc;

%% INITIALIZATION
%% Simulation parameters
% prepare sensor data to fed to the net
run_steps          = 1500;
% iterations to present a data point to the net
show_step          = 1;
% runtime parameters
run_iters          = run_steps - 1;
run_iters_extended = run_iters;
% number of maps
maps               = 2;

% id of the current update rule
rand_update_rule  = 1;
% possible update rules - depends on net topology
update_rules      = 2;
% individual id for each rule
rules_ids         = 1:update_rules;
% convergence steps counter
convergence_steps = 1;
% simulation points (sensor data is presented for half simulation)
sim_points        = run_iters_extended;
%% Network elements setup

% network iterator
net_iter = 1;

% init maps and indices
m1 = rand; m2 = rand;
% maps ids
m1_id = 1; m2_id = 2;

% temporal elements
m1_ant = m1; m1_ant_ant = m1_ant; m2_ant = m2;
int_m1      = 0.0;
deriv_m1    = 0.0;
deriv2_m1   = 0.0;
deriv_m2    = 0.0;
dtk         = 0.01;

% number of connections of the maps (sensor + relations to other maps)
m1_links = 1:1; m2_links = 1:1;

% init maps errors for learning rate adaptation
em1 = zeros(1, length(m1_links)); em2 = zeros(1, length(m2_links));


% learning rates setup
maps_nr = maps;
sensors_nr = maps_nr;
error_nr = update_rules;
lrates_nr = update_rules;
net_data = zeros(sim_points, maps_nr+error_nr+lrates_nr);

% init learning rates and bounds
ETA = 0.002;
ETAH = 10*ETA;

etam1 = ETA*ones(sim_points, length(m1_links));
etam2 = ETA*ones(sim_points, length(m2_links));


%% NETWORK DYNAMICS
%
% encoded relationships
%
%   m2 = int(m1) |--> m1 = dm2/dt

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
    
    %% go through current convergence step update sequence
    
    for i = 1:update_rules
        % pick a rule to update
        rand_update_rule = rules_ids(i);
        % update the corresponding node using simple error
        switch rand_update_rule
            % Circular permutations of the error signal
            %% m1 update
            case 1
                m1 = (1-etam1(convergence_steps, m1_links(1)))*m1 + etam1(convergence_steps, m1_links(1)) * ...
                    deriv_m2;
                %% m2 update
            case 2
                m2 = (1-etam2(convergence_steps, m2_links(1)))*m2 + etam2(convergence_steps, m2_links(1)) * ...
                    int_m1;
                
        end
        
        % compute the temporal representations
        deriv2_m1   = (m1 - 2*m1_ant + m1_ant_ant)/(dtk^2);
        int_m1      = int_m1 + (m1_ant + m1)*dtk*0.5;
        deriv_m2    = (m2 - m2_ant)/dtk;
        
        
    end
    % --------errors--------
    % error computation for learning rate adaptation
    % m1
    em1(1) = m1 - deriv_m2;
    % m2
    em2(1) = m2 - int_m1;
    
    
    %% WRITE DATA TO STRUCT
    % maps
    net_data(convergence_steps, 1) = m1;
    net_data(convergence_steps, 2) = m2;
    
    % error
    net_data(convergence_steps, 3) = em1(1);
    net_data(convergence_steps, 4) = em2(1);
    
    % sample idx
    net_data(convergence_steps, 5) = convergence_steps;
    
    % add the learning rates in the struct
    net_data(convergence_steps, 6) = etam1(convergence_steps, m1_links(1));
    net_data(convergence_steps, 7) = etam2(convergence_steps, m2_links(1));
    
    
    %% update history and loop indices
    m1_ant_ant = m1_ant; m1_ant = m1; m2_ant = m2;
    net_iter = net_iter + 1;
    convergence_steps = convergence_steps + 1;
end

% fill in the net simulation data into the visualization struct
fusion_analyzer_data = net_data;

%% VISUALIZATION
figure(1);
% ---------------- R1 ----------------------

plot(net_data(:, 1), '.b');hold on;
plot(net_data(:, 2), '.r');
legend('m_i', 'm_j');
xlabel('iterations'); ylabel('value');
% -------------------------------------------
% set figure props
set(gcf,'color','w');
box off;
% -------------- Erros signals --------------
% figure(2);
% %-------------------M1-------------------
% he(1) = subplot(2,1,1);
% plot(net_data(:,3),'.b');
% title('m_i error');
% grid on;
% %-------------------M2-------------------
% he(2) = subplot(2,1,2);
% plot(net_data(:,4), '.r');
% title('m_j error');grid on;
% % set figure props
% set(gcf,'color','w');

