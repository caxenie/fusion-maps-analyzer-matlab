function l = update_learning_rate(oldl, err, idx, ~, type)
switch type
    case 'divisive'
        sum = 0;
        default_val = 0.02;
        for k = 1:(length(err))
            if(k~=idx)
                sum = sum + abs(err(k));
            end
        end
        l = default_val*(sum/((length(err)-1)*abs(err(idx))));
    case 'decay'
        sum = 0;
        default_val = 0.02;
        for k = 1:(length(err))
            if(k~=idx)
                sum = sum + abs(err(k));
            end
        end
        l = default_val*(sum/((length(err)-1)*abs(err(idx))));
        l = oldl + l;   
end
end