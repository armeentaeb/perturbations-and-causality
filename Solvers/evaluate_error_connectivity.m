function [xaxis,yaxis] = evaluate_error_connectivity(B)
global Bstar
p = size(Bstar,1);
Bvec = reshape(B,p^2,1);
Bstar_vec = reshape(Bstar,p^2,1);

[sort_val, sort_ind] = sort(abs(Bvec),'descend');
tot_nonzeros = length(find(abs(Bvec)>=10^(-3)));
xaxis = zeros(p^2,1);
yaxis = zeros(p^2,1);
agg_y = 0;
agg_x = 0;

for i = 1:tot_nonzeros
    if abs(Bstar_vec(sort_ind(i))) >= 10^(-3)
        xaxis(i) = agg_x;
        if i == 1
            yaxis(i) = 1;
            agg_y = 1;
        else
            yaxis(i) = agg_y+1;
            agg_y = agg_y+1;
        end
    else
        yaxis(i) = agg_y;
        if i == 1
            xaxis(i) = 1;
            agg_x = 1;

        else
            xaxis(i) = agg_x+1;
            agg_x = agg_x+1;
        end
    end
    
end
xaxis(i+1:end) = xaxis(i);
yaxis(i+1:end) = yaxis(i);


end