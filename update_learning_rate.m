function l = update_learning_rate(oldl, err, idx, default_val, type)
switch type
    case 'fixed'
        l = default_val;
    case 'divisive'
        sum = 0;
        for k = 1:(length(err))
            if(k~=idx)
                sum = sum + abs(err(k));
            end
        end
        l = (sum/((length(err)-1)*abs(err(idx))));
    case 'decay'
        sum = 0;
        for k = 1:(length(err))
            if(k~=idx)
                sum = sum + abs(err(k));
            end
        end
        l = default_val*(sum/((length(err)-1)*abs(err(idx))));
        l = oldl + l;
end
end