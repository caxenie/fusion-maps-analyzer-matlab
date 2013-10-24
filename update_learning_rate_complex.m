function l = update_learning_rate_complex(l_old, grad_bar_old, derr, derr_old, u, d, l_min, l_max, k, gama, type)
switch type
    case 'up-down-factor'
        if(derr*derr_old > 0)
            l = min(l_old*u, l_max);
        else if(derr*derr_old < 0)
                l = max(l_old*d, l_min);
            else
                l = l_old;
            end
        end
    case 'delta-bar-delta'
        if(grad_bar_old*derr > 0)
            l = l_old + k;
        else if(grad_bar_old*derr < 0)
                l = (1-gama)*l_old;
            else
                l = l_old;
            end
        end
end
end