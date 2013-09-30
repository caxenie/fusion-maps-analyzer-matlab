% Fusion network for testing simple algebraic and temporal relationships
close all;
clear all;
clc;

% load data 
input_samples = 5;
maps = 6;                       % number of maps
input_data = rand([input_samples, maps]);
m1_sensor = input_data(:, 1);
m2_sensor = input_data(:, 2);
m3_sensor = input_data(:, 3);
m4_sensor = input_data(:, 4);
m5_sensor = input_data(:, 5);
m6_sensor = input_data(:, 6);

%% TODO generate some datasets with data in [-1 1]


%% Initialization 
rand_update_rule = 1;           % id of the current update rule
update_rules = 14;              % possible update rules - depends on net topology
rules_ids = 1:update_rules;     % individual id for each rule
convergence_steps = 1;          % convergence steps counter
data_show_freq = 10000;           % freq to show data point to net
sim_points = length(input_data)*data_show_freq;

% sensor connection state - connected to net / not
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

% net params
% learning rates setup 
maps_nr = maps;
error_nr = update_rules;
lrates_nr = update_rules;
net_data = zeros(sim_points, maps_nr+error_nr+lrates_nr);

% init learning rates 
ETA = 0.002;
ETAH = 0.02;

etam1 = ETA*ones(sim_points, length(m1_links)); 
etam2 = ETA*ones(sim_points, length(m2_links));
etam3 = ETA*ones(sim_points, length(m3_links)); 
etam4 = ETA*ones(sim_points, length(m4_links));
etam5 = ETA*ones(sim_points, length(m5_links)); 
etam6 = ETA*ones(sim_points, length(m6_links));

%% Network dynamics
% encoded relationships
%   m2 = 3*m1
%   m3 = m2*m4
%   m4 = m5 + 2*m6

while(1)
    % check for end of simulation 
    if(net_iter == length(input_data))
        if (net_iter*data_show_freq == convergence_steps)
            break;
        end;
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
       em4(3) = m4 - (m5+2*m6);
       % m5
       em5(1) = m5 - m5_sensor(net_iter);
       em5(2) = m5 - (m4 - 2*m6);
       % m6 
       em6(1) = m6 - m6_sensor(net_iter);
       em6(2) = m6 - ((m4 - m5)/2); 
    end
    
    %% TODO update learing rates 
    
    % write data to net struct 
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
    
    % learning rates - fixed for the moment
    
    net_data(convergence_steps, 21) = convergence_steps;
    
    % check new sensor data presentation
    if(convergence_steps == (net_iter+1)*data_show_freq - 1)
        net_iter = net_iter + 1;
    end
    convergence_steps = convergence_steps + 1;
end

%% Visualization 
figure(1);
subplot(3,1,1);
plot(net_data(:,1), '.r');
subplot(3,1,2);
plot(net_data(:,2), '.b');
subplot(3,1,3);
plot(net_data(:,1), net_data(:,2));
