function y = clamp(a,b)
    if(a - b > 0.00000001)
        y = b;
    else
        y = a;
    end
end