function l = update_confidence_factor(oldl, err, idx, default_val, type)
switch type
    case 'incdec'
        u = 0.00000001;
        sum = 0;
        for k = 1:(length(err))
            if(k~=idx)
                sum = sum + abs(err(k));
            end
        end
        if((sum/((length(err)-1)*abs(err(idx))))>1)
            l = oldl + u;
        else
            if((sum/((length(err)-1)*abs(err(idx))))<1)
                l = oldl - u;
            else
                l = oldl;
            end
        end
    case 'divisive'
        sum = 0;
        for k = 1:(length(err))
            if(k~=idx)
                sum = sum + abs(err(k));
            end
        end
        l = default_val * (sum/((length(err)-1)*abs(err(idx))));
    case 'decay'
        sum = 0;
        for k = 1:(length(err)) 
            if(k~=idx)
                sum = sum + abs(err(k));
            end
        end
        l = oldl + default_val*(sum/((length(err)-1)*abs(err(idx))));   
    case 'fixed'
        l = default_val;
end
end