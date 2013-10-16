function l = update_learning_rate(oldl, err, idx, default_val, type)
switch type
    case 'divisive'
        sum = 0;
        for k = 1:(length(err))
            if(k~=idx)
                sum = sum + abs(err(k));
            end
        end
        l = default_val*(((length(err)-1)*abs(err(idx)))/sum);
    case 'decay'
        sum = 0;
        for k = 1:(length(err))
            if(k~=idx)
                sum = sum + abs(err(k));
            end
        end
        l = default_val*(((length(err)-1)*abs(err(idx)))/sum);
        l = oldl + l;
end
end