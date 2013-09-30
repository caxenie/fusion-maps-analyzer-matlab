function l = update_learning_rate(err, r, maps)
    sum = 0;
    ETA = 0.02;
        for i = 1:maps
            if(i~=r)
                sum = sum + abs(err(i));
            end
        end
    l = ETA*(sum/(maps-1)*abs(err(r)));
end