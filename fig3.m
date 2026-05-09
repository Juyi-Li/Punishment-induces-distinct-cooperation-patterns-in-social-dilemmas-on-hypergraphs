e = [];
for kk = 0.01:0.01:0.99
    for jj = 0.01:0.01:0.99
        if 1-kk-jj>=0.01
            e = [e;kk,jj,1-kk-jj];
        end
    end
end

% Parameter definitions
R = 1; S = -0.1; T = 1.8; P = 0; G = 0.8; W = 0.1; theta = 0.7;
tmax = 4000;
dt = 0.1;

% Generate all initial conditions using the precomputed simplex grid e
N0 = size(e,1);          % Number of initial conditions
init = e;                % N0 × 3 matrix, each row is an initial condition [y1, y2, y3]

% Initialize the result array
CE_domain = zeros(101, 101, 101);

% Vectorized and parallelized loops over alpha, beta, and delta
parfor alpha_idx = 1:101
    alpha = (alpha_idx - 1) * 0.01;
    
    % Preallocate local results for each alpha to avoid write conflicts in parfor
    local_CE = zeros(101, 101);
    
    for beta_idx = 1:101
        beta = (beta_idx - 1) * 0.01;
        
        for delta_idx = 1:101
            tic
            delta = (delta_idx - 1) * 0.01;
            
            % Current state matrix: N0 × 3
            Y = init;   % Each row is [y1, y2, y3]
            
            % Time evolution, retaining only the final state
            for t = 1:tmax
                y1 = Y(:,1);   % N0 × 1 vector
                y2 = Y(:,2);
                y3 = 1 - y1 - y2;
                
                % Compute pi_c, pi_d, and pi_p in a vectorized form
                % Note: all operations are performed on column vectors
                pi_c = (1-delta)*( y1*(R/2*(1+theta)) + y2*S + y3*(R/2*(1+theta)) ) + ...
                       delta*( y1.^2*(R/3*(1+theta+theta^2)) + 2*y1.*y2*(G/2*(1+theta)) + 2*y1.*y3*(R/3*(1+theta+theta^2)) + ...
                               y2.^2*S + 2*y2.*y3*(G/2*(1+theta)) + y3.^2*(R/3*(1+theta+theta^2)) );
                
                pi_d = (1-delta)*( y1*T + y2*P + y3*(T-beta) ) + ...
                       delta*( y1.^2*(T/2*(1+theta)) + 2*y1.*y2*W + 2*y1.*y3*(T/2*(1+theta)-beta) + ...
                               y2.^2*P + 2*y2.*y3*(W-beta) + y3.^2*(T/2*(1+theta)-2*beta) );
                
                pi_p = (1-delta)*( y1*(R/2*(1+theta)) + y2*(S-alpha) + y3*(R/2*(1+theta)) ) + ...
                       delta*( y1.^2*(R/3*(1+theta+theta^2)) + 2*y1.*y2*(G/2*(1+theta)-alpha) + 2*y1.*y3*(R/3*(1+theta+theta^2)) + ...
                               y2.^2*(S-2*alpha) + 2*y2.*y3*(G/2*(1+theta)-alpha) + y3.^2*(R/3*(1+theta+theta^2)) );
                
                pi_bar = y1.*pi_c + y2.*pi_d + y3.*pi_p;
                
                % Update y1 and y2
                Y(:,1) = y1 + dt * y1 .* (pi_c - pi_bar);
                Y(:,2) = y2 + dt * y2 .* (pi_d - pi_bar);
                
                % Boundary handling
                Y(:,1) = max(0, min(1, Y(:,1)));
                Y(:,2) = max(0, min(1, Y(:,2)));
            end  % End of time evolution
            
            % Compute y1 + y3 for each initial condition at the final state
            y3_final = 1 - Y(:,1) - Y(:,2);
            result = mean( Y(:,1) + y3_final );   % Average over all initial conditions
            
            local_CE(beta_idx, delta_idx) = result;
            toc
        end
    end
    
    % Write the local results into the global array
    % The indexing in parfor must be simple
    CE_domain(alpha_idx, :, :) = local_CE;
end