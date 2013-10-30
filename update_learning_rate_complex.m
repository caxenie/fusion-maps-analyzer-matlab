function l = update_learning_rate_complex(l_old, grad_bar_old, derr, derr_old, u, d, l_min, l_max, k, gama, type)
switch type
    case 'grad-history'
        if(derr*derr_old > 0)
            l = min(l_old*u, l_max);
        end
        if(derr*derr_old < 0)
                l = max(l_old*d, l_min);
        end
        if(derr*derr_old==0)
            l = l_old;
        end
    case 'grad-history-avergae'
        if(grad_bar_old*derr > 0)
            l = l_old + k;
        end
        if(grad_bar_old*derr < 0)
                l = (1-gama)*l_old;
        end
        if(grad_bar_old*derr ==0)
                l = l_old;
        end
end
end