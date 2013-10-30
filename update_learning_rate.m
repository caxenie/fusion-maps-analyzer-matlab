function l = update_learning_rate(oldl, err, idx, default_val, type)
switch type
    case 'fixed'
        l = default_val;
    case 'adaptive'
        u = 1.5;
        sum = 0;
        for k = 1:(length(err))
            if(k~=idx)
                sum = sum + (err(k));
            end
        end
        if((sum/((length(err)-1)*(err(idx))))>1)
            l = oldl * u;
        end
        if((sum/((length(err)-1)*(err(idx))))<1)
            l = oldl * 1/u;
        end
        if((sum/((length(err)-1)*(err(idx))))==1)
            l = oldl;
        end
    case 'divisive'
        sum = 0;
        for k = 1:(length(err))
            if(k~=idx)
                sum = sum + (err(k));
            end
        end
        l = default_val*(sum/((length(err)-1)*(err(idx))));   
    case 'decay'
        sum = 0;
        for k = 1:(length(err))
            if(k~=idx)
                sum = sum + (err(k));
            end
        end
        l = default_val*(sum/((length(err)-1)*(err(idx))));
        l = oldl + l;
end
end