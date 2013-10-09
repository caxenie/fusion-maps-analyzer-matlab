function l = update_learning_rate_complex(l_old, derr, derr_old, u, d, l_min, l_max, type)
switch type
    case 'up-down-factor'
        if(derr/derr_old>=1)
            l = l_old*u;
        else
            l = l_old*d;
        end
    case 'min-max-factor'
        if(derr/derr_old > 1)
            l = min(l_old*u, l_max);
        else if(derr/derr_old < 1)
                l = max(l_old*d, l_min);
            else
                l = l_old;
            end
        end
end
end