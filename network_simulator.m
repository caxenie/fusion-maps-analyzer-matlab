%% Fusion network for testing simple algebraic and temporal relationships
close all;
clear all;
clc;

% prepare sensor data  to fed to the net
run_iters = 320000;
maps      = 6;                       % number of maps
m1_sensor = rand;
m2_sensor = rand;
m3_sensor = rand;
m4_sensor = rand;
m5_sensor = rand;
m6_sensor = rand;



%% Initialization 
rand_update_rule = 1;           % id of the current update rule
update_rules = 14;              % possible update rules - depends on net topology
rules_ids = 1:update_rules;     % individual id for each rule
convergence_steps = 1;          % convergence steps counter
sim_points = run_iters;
% clamp connection time
clamp1_on = 20000;
clamp2_on = 8000;
clamp1_off = 300000;
clamp2_off = 18000;
clamp_value = 1;

% sensor connection state - connected to net / not
sensor_connected = zeros(1, maps);

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
net_data = zeros(sim_points, maps_nr+error_nr);

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
    if (convergence_steps == run_iters)
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
    
    net_data(convergence_steps, 21) = convergence_steps;
       
    % TODO learning rates - fixed for the moment
 
    % sensor values clamping 
    switch(convergence_steps)
        case clamp1_on
            sensor_connected(m1_id)  = 1;
            m1_sensor = clamp_value;
        case clamp1_off
            sensor_connected(m1_id)  = 0;
        case clamp2_on
            sensor_connected(m2_id)  = 1;
            m2_sensor = -clamp_value;
        case clamp2_off
            sensor_connected(m2_id) = 0;
    end
    
    convergence_steps = convergence_steps + 1;
end

%% Visualization 
figure(1);
% ---------------- R1 ----------------------
subplot(4, 3, 1);
plot(net_data(:, 21), net_data(:, 1), '.r');
grid on; title('M1 map values');

subplot(4, 3, 4);
plot(net_data(:, 21), net_data(:, 2), '.b');
grid on; title('M2 map values');

subplot(4, 3, 10);
plot(net_data(:, 2), net_data(:, 1), '.k');
grid on; title('M1 map Vs. M2 map values');
% -------------------------------------------

% ---------------- R2 ----------------------
subplot(4, 3, 2);
plot(net_data(:, 21), net_data(:, 2), '.r');
grid on; title('M2 map values');

subplot(4, 3, 5);
plot(net_data(:, 21), net_data(:, 3), '.b');
grid on; title('M3 map values');

subplot(4, 3, 8);
plot(net_data(:, 21), net_data(:, 4), '.g');
grid on; title('M4 map values');

subplot(4, 3, 11);
plot3(net_data(:, 2), net_data(:, 3), net_data(:, 4));
grid on; title('M2 map Vs. M3 map values Vs. M4 values');
% -------------------------------------------

% ---------------- R2 ----------------------
subplot(4, 3, 3);
plot(net_data(:, 21), net_data(:, 4), '.r');
grid on; title('M2 map values');

subplot(4, 3, 6);
plot(net_data(:, 21), net_data(:, 5), '.b');
grid on; title('M3 map values');

subplot(4, 3, 9);
plot(net_data(:, 21), net_data(:, 6), '.g');
grid on; title('M4 map values');

subplot(4, 3, 12);
plot3(net_data(:, 4), net_data(:, 5), net_data(:, 6));
grid on; title('M4 map Vs. M5 map values Vs. M6 values');
% -------------------------------------------


% % fill in the data 
% fusion_analyzer_data = net_data;
% %--------------------------------------------------------------------------------------------------------
% fig2Handle = figure(2);
% set(fig2Handle, 'Position', [100, 100, 600, 1000]);
% %-------------------------
% subplot(3,1,1)
% % axis control and adjustment
% [ax1] = axescheck(fusion_analyzer_data(:,20), fusion_analyzer_data(:,1));
% ax1 = newplot(ax1);
% set(ax1,'XGrid','on');
% set(ax1,'YGrid','on');
% if ~ishold(ax1)
%   [minx1,maxx1] = minmax(fusion_analyzer_data(:,20));
%   [miny1,maxy1] = minmax(fusion_analyzer_data(:,1));
%   axis(ax1,[minx1 maxx1 miny1 maxy1])
% end
% title('Maps values during simulation');
% ylabel('Map 1 values')
% xlabel('Data points')
% hc1 =  line('parent',ax1, 'linestyle','-','erase','xor','xdata',fusion_analyzer_data(1,20),'ydata',fusion_analyzer_data(1,1));
% ha1 = line('parent',ax1,'linestyle','-','erase','none','xdata',[],'ydata',[]);
% %-------------------------
% subplot(3,1,2)
% % axis control and adjustment
% [ax2] = axescheck(fusion_analyzer_data(:,20), fusion_analyzer_data(:,2));
% ax2 = newplot(ax2);
% set(ax2,'XGrid','on');
% set(ax2,'YGrid','on');
% if ~ishold(ax2)
%   [minx2,maxx2] = minmax(fusion_analyzer_data(:,20));
%   [miny2,maxy2] = minmax(fusion_analyzer_data(:,2));
%   axis(ax2,[minx2 maxx2 miny2 maxy2])
% end
% ylabel('Map 2 values')
% xlabel('Data points')
% hc2 =  line('parent',ax2,'linestyle','-','erase','xor','xdata',fusion_analyzer_data(1,20),'ydata',fusion_analyzer_data(1,2));
% ha2 = line('parent',ax2,'linestyle','-','erase','none','xdata',[],'ydata',[]);
% %-------------------------
% subplot(3,1,3)
% % axis control and adjustment
% [ax3] = axescheck(fusion_analyzer_data(:,1), fusion_analyzer_data(:,2));
% ax3 = newplot(ax3);
% if ~ishold(ax3)
%   [minx3,maxx3] = minmax(fusion_analyzer_data(:,1));
%   [miny3,maxy3] = minmax(fusion_analyzer_data(:,2));
%   axis(ax3,[minx3 maxx3 miny3 maxy3])
% end
% t = -1.0:1.0;
% plot(t, 3*t, '-m');
% set(ax3,'XGrid','on');
% set(ax3,'YGrid','on');
% hc3 =  line('parent',ax3,'linestyle','-','erase','xor','xdata',fusion_analyzer_data(1,1),'ydata',fusion_analyzer_data(1,2));
% ha3 = line('parent',ax3,'linestyle','-','erase','none','xdata',[],'ydata',[]);
% title('M2 dependency on M1');
% xlabel('M1 values');
% ylabel('M2 values');
% % small pause to start recording with external device
% pause(5); % s
% % dynamic plot simulataneously in all subplots
% for i=2:length(fusion_analyzer_data(:,20))
%     j = i-1:i;
%     set(hc1,'xdata', fusion_analyzer_data(i,20), 'ydata', fusion_analyzer_data(i,1));
%     set(ha1,'xdata', fusion_analyzer_data(j,20), 'ydata', fusion_analyzer_data(j,1));
%     set(hc2,'xdata', fusion_analyzer_data(i,20), 'ydata', fusion_analyzer_data(i,2));
%     set(ha2,'xdata', fusion_analyzer_data(j,20), 'ydata', fusion_analyzer_data(j,2));
%     set(hc3,'xdata', fusion_analyzer_data(i,1), 'ydata', fusion_analyzer_data(i,2));
%     set(ha3,'xdata', fusion_analyzer_data(j,1), 'ydata', fusion_analyzer_data(j,2));
%     drawnow;
% end
% 
% %--------------------------------------------------------------------------------------------------------
% fig3Handle = figure(3);
% set(fig3Handle, 'Position', [900, 900, 600, 1200]);
% subplot(4,1,1)
% % axis control and adjustment
% [ax4] = axescheck(fusion_analyzer_data(:,20), fusion_analyzer_data(:,3));
% ax4 = newplot(ax4);
% set(ax4,'XGrid','on');
% set(ax4,'YGrid','on');
% if ~ishold(ax4)
%   [minx4,maxx4] = minmax(fusion_analyzer_data(:,20));
%   [miny4,maxy4] = minmax(fusion_analyzer_data(:,3));
%   axis(ax4,[minx4 maxx4 miny4 maxy4])
% end
% title('Maps values during simulation');
% hc4 =  line('parent',ax4, 'linestyle','-','erase','xor','xdata',fusion_analyzer_data(1,20),'ydata',fusion_analyzer_data(1,3));
% ha4 = line('parent',ax4,'linestyle','-','erase','none','xdata',[],'ydata',[]);
% ylabel('Map 3 values')
% xlabel('Data points')
% 
% %-------------
% subplot(4,1,2)
% % axis control and adjustment
% [ax5] = axescheck(fusion_analyzer_data(:,20), fusion_analyzer_data(:,4));
% ax5 = newplot(ax5);
% set(ax5,'XGrid','on');
% set(ax5,'YGrid','on');
% if ~ishold(ax5)
%   [minx5,maxx5] = minmax(fusion_analyzer_data(:,20));
%   [miny5,maxy5] = minmax(fusion_analyzer_data(:,4));
%   axis(ax5,[minx5 maxx5 miny5 maxy5])
% end
% hc5 =  line('parent',ax5, 'linestyle','-','erase','xor','xdata',fusion_analyzer_data(1,20),'ydata',fusion_analyzer_data(1,4));
% ha5 = line('parent',ax5,'linestyle','-','erase','none','xdata',[],'ydata',[]);
% ylabel('Map 4 values')
% xlabel('Data points')
% %-------------
% subplot(4,1,3)
% % axis control and adjustment
% [ax6] = axescheck(fusion_analyzer_data(:,20), fusion_analyzer_data(:,5));
% ax6 = newplot(ax6);
% set(ax6,'XGrid','on');
% set(ax6,'YGrid','on');
% if ~ishold(ax6)
%   [minx6,maxx6] = minmax(fusion_analyzer_data(:,20));
%   [miny6,maxy6] = minmax(fusion_analyzer_data(:,5));
%   axis(ax6,[minx6 maxx6 miny6 maxy6])
% end
% hc6 = line('parent',ax6, 'linestyle','-','erase','xor','xdata',fusion_analyzer_data(1,20),'ydata',fusion_analyzer_data(1,5));
% ha6 = line('parent',ax6,'linestyle','-','erase','none','xdata',[],'ydata',[]);
% ylabel('Map 5 values')
% xlabel('Data points')
% %------------
% subplot(4,1,4)
% % axis control and adjustment
% [ax7] = axescheck(fusion_analyzer_data(:,20), fusion_analyzer_data(:,6));
% ax7 = newplot(ax7);
% set(ax7,'XGrid','on');
% set(ax7,'YGrid','on');
% if ~ishold(ax7)
%   [minx7,maxx7] = minmax(fusion_analyzer_data(:,20));
%   [miny7,maxy7] = minmax(fusion_analyzer_data(:,6));
%   axis(ax7,[minx7 maxx7 miny7 maxy7])
% end
% hc7 =  line('parent',ax7, 'linestyle','-','erase','xor','xdata',fusion_analyzer_data(1,20),'ydata',fusion_analyzer_data(1,6));
% ha7 = line('parent',ax7,'linestyle','-','erase','none','xdata',[],'ydata',[]);
% ylabel('Map 6 values')
% xlabel('Data points')
% % small pause to start recording with external device
% pause(5); % s
% % loop and dynamically shouw map convergence
% for i=2:length(fusion_analyzer_data(:,20))
%     j = i-1:i;
%     set(hc4,'xdata', fusion_analyzer_data(i,20), 'ydata', fusion_analyzer_data(i,3));
%     set(ha4,'xdata', fusion_analyzer_data(j,20), 'ydata', fusion_analyzer_data(j,3));
%     set(hc5,'xdata', fusion_analyzer_data(i,20), 'ydata', fusion_analyzer_data(i,4));
%     set(ha5,'xdata', fusion_analyzer_data(j,20), 'ydata', fusion_analyzer_data(j,4));
%     set(hc6,'xdata', fusion_analyzer_data(i,20), 'ydata', fusion_analyzer_data(i,5));
%     set(ha6,'xdata', fusion_analyzer_data(j,20), 'ydata', fusion_analyzer_data(j,5));
%     set(hc7,'xdata', fusion_analyzer_data(i,20), 'ydata', fusion_analyzer_data(i,6));
%     set(ha7,'xdata', fusion_analyzer_data(j,20), 'ydata', fusion_analyzer_data(j,6));
%     drawnow;
% end
% 
% %--------------------------------------------------------------------------------------------------------
% % per relation errors
% % first relation between M1 and M2
% figure(3);
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
% figure(4);
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
% figure(5);
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
% %--------------------------------------------------------------------------------------------------------
% % analize the first relationship R1: M2 = 3*M1
% % analize the (M1, M2) dependency 
% figure(6);
% plot(fusion_analyzer_data(:,20), 3*fusion_analyzer_data(:,1));
% hold on;
% plot(fusion_analyzer_data(:,20), fusion_analyzer_data(:,2), '-r');
% title('The dependency between M1 and M2');
% grid on;
% legend('3*M1 data', 'M2 data');
% xlabel('Data points');
% 
% %--------------------------------------------------------------------------------------------------------
% % analize the second relationship R2: M3 = M2*M4
% % analize the (M2, M3, M4) dependency 
% figure(7);
% plot(fusion_analyzer_data(:,20), fusion_analyzer_data(:,2).*fusion_analyzer_data(:,4));
% hold on;
% plot(fusion_analyzer_data(:,20), fusion_analyzer_data(:,3), '-r');
% title('The dependency between M2, M3, and M4');
% legend('M2*M4 data', 'M3 data');
% grid on;
% xlabel('Data points');
% 
% %--------------------------------------------------------------------------------------------------------
% % analize the second relationship R3: M4 = M5 + 2*M6
% % analize the (M4, M5, M6) dependency 
% figure(8);
% plot(fusion_analyzer_data(:,20), fusion_analyzer_data(:,5)+(2*fusion_analyzer_data(:,6)));
% hold on;
% plot(fusion_analyzer_data(:,20), fusion_analyzer_data(:,4), '-r');
% title('The dependency between M4, M5, and M6');
% legend('M5+2*M6 data', 'M4 data');
% grid on;
% xlabel('Data points');  
