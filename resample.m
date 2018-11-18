function [xu_out, q_out] = resample(xu_in, q_in)
    [nl,nc] = size(xu_in);
    N = nc;
    edges_vec = [0 cumsum(q_in)];       % Distribution of cumulative weights
    rand_shoot = rand(N,1);             % Uniform distribution
    sorted_shoot = sort(rand_shoot);    % Sorting of the results
    % Counts how many randomly shot points fall into the cumulative distribution
    % The more point fall to a bin, the more likely the real measure is there
    [count, bins] = histc(sorted_shoot,edges_vec);     
    xu_out = xu_in(:,bins);             % Populates the bins accordingly
    q_out = ones(1,N)/N;                % Gives uniform weight to particles
end