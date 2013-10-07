% build a sensor input data set
function y = build_data_set(amplitude_step, switch_time, size)
    y(1:switch_time) = rand;
    for i=switch_time+1:size
        if(mod(i,switch_time)==0)
            amplitude_step = 2*amplitude_step; 
        end
        y(i) = amplitude_step;
    end
end