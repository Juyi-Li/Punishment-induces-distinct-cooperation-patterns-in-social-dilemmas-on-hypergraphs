%% ====== CDP: rewritten according to the new CD update process (three strategies C/D/P) ======
% Description:
% 1) At each step, a focal node i is randomly selected.
% 2) If node i has both pairwise and three-body neighbors, the update channel is selected with equal probability.
% 3) If only one type of neighbor exists, that channel is used.
% 4) Payoffs are still computed using the original CDP pair_avg and triple_avg terms, then weighted by delta.
% 5) The focal node i updates its strategy according to the Fermi rule.
%
% State encoding:
%   D = 0
%   C = 2
%   P = 7

clearvars -except A2 H3 ER_lower_order ER_higher_order

%% -------------------- Choose input matrices --------------------
if exist('A2','var') && exist('H3','var')
    B = A2;
    M = H3;
elseif exist('ER_lower_order','var') && exist('ER_higher_order','var')
    B = ER_lower_order;
    M = ER_higher_order;
else
    error('Please provide A2 and H3, or ER_lower_order and ER_higher_order in the workspace first.');
end

rng(42);

%% -------------------- Check and clean the network --------------------
B = spones(B);
n = size(B,1);

B = triu(B,1) + triu(B,1)';              % Force an undirected structure without self-loops
assert(isequal(B,B.'), 'B must be symmetric, corresponding to an undirected graph.');
assert(all(diag(B)==0), 'B must not contain self-loops; the diagonal entries should be zero.');

M = spones(M);
if size(M,1) ~= n
    error('The number of rows in M (%d) must be consistent with the size of B (%d).', size(M,1), n);
end

colsum = full(sum(M,1));
keep   = (colsum == 3);
if any(~keep)
    warning('Removed %d non-3-uniform hyperedge columns.', nnz(~keep));
    M = M(:, keep);
end
H = size(M,2);

fprintf('n = %d\n', n);
fprintf('2-body edges = %d\n', nnz(triu(B,1)));
fprintf('3-body hyperedges = %d\n', H);

%% -------------------- Precompute pairwise neighborhoods and three-body hyperedges --------------------
Nhood = cell(n,1);
for i = 1:n
    Nhood{i} = find(B(i,:));
end

edgeNodes = cell(H,1);
for h = 1:H
    v = find(M(:,h))';
    if numel(v) == 3
        edgeNodes{h} = v;
    else
        edgeNodes{h} = [];
    end
end

nodeEdges = cell(n,1);
for i = 1:n
    nodeEdges{i} = find(M(i,:));
end

deg2 = cellfun(@numel, Nhood);
deg3 = cellfun(@numel, nodeEdges);

fprintf('nodes with 2-body neighbors: %d / %d\n', nnz(deg2>0), n);
fprintf('nodes with 3-body neighbors: %d / %d\n', nnz(deg3>0), n);

%% ===================== Dynamics and payoff parameters =====================
omega = 0.1;

S     = -0.1;
T     = 1.8;
G     = 0.8;
W     = 0.1;
alpha = 0.7;
beta  = 0.7;

% Delta still represents the weight of the three-body payoff contribution
deltas = 0.0:0.1:1.0;

rep    = 100;
tmax   = 500000;

checkEvery = 500;
bufLen     = 10;
tol_std    = 1e-4;

%% ===================== Construct the initial composition grid for (C,D,P), excluding zero components =====================
step = 0.1;
comps = [];

for c = 0:step:1
    for d = 0:step:(1-c)
        p = 1 - c - d;
        if p < -1e-9
            continue;
        end
        if abs(round(c/step) - c/step) > 1e-9 || ...
           abs(round(d/step) - d/step) > 1e-9 || ...
           abs(round(p/step) - p/step) > 1e-9
            continue;
        end
        if c <= 0 || d <= 0 || p <= 0
            continue;
        end
        comps = [comps; c, d, p]; %#ok<AGROW>
    end
end

nComps = size(comps,1);
fprintf('Total initial (C,D,P) compositions: %d\n', nComps);

%% ===================== Result containers =====================
Dnum   = numel(deltas);

C_end  = zeros(rep, nComps, Dnum);
P_end  = zeros(rep, nComps, Dnum);
CP_end = zeros(rep, nComps, Dnum);

meanC  = zeros(nComps, Dnum);
meanP  = zeros(nComps, Dnum);
meanCP = zeros(nComps, Dnum);

%% ===================== Main loops =====================
for dIdx = 1:Dnum
    delta = deltas(dIdx);
    fprintf('\n===== delta = %.2f =====\n', delta);

    for compIdx = 1:nComps
        cFrac = comps(compIdx,1);
        dFrac = comps(compIdx,2);
        pFrac = comps(compIdx,3);

        nC = round(cFrac * n);
        nD = round(dFrac * n);
        nP = n - nC - nD;    % Use P to absorb rounding errors

        parCend = zeros(rep,1);
        parPend = zeros(rep,1);

        parfor r = 1:rep
            tic
            % ---------- Initialization ----------
            state = zeros(n,1,'int8');   % D = 0

            if nC > 0
                idxC = randperm(n, nC);
                state(idxC) = 2;         % C
            end

            if nP > 0
                remain = find(state==0);
                nP_loc = min(nP, numel(remain));
                if nP_loc > 0
                    idxP = remain(randperm(numel(remain), nP_loc));
                    state(idxP) = 7;     % P
                end
            end

            % Rolling buffer for early stopping
            cbuf = nan(bufLen,1);
            ptr  = 0;

            for t = 1:tmax
                old = state;

                % 1) Pick a focal node i uniformly at random
                i = randi(n);

                nei_i = Nhood{i};
                Hi    = nodeEdges{i};

                has2_i = ~isempty(nei_i);
                has3_i = ~isempty(Hi);

                % If i has neither pairwise nor three-body neighbors, skip this step
                if ~(has2_i || has3_i)
                    continue;
                end

                % 2) Choose the pairwise or three-body update channel with equal probability
                if has2_i && has3_i
                    use3 = (rand < 0.5);   % 50%-50%
                elseif has2_i
                    use3 = false;
                else
                    use3 = true;
                end

                % 3) Choose the role model j
                if ~use3
                    % Choose j from pairwise neighbors
                    j = nei_i(randi(numel(nei_i)));
                else
                    % First choose one incident 3-hyperedge
                    h = Hi(randi(numel(Hi)));
                    v = edgeNodes{h};      % 1 x 3
                    others = v(v ~= i);
                    j = others(randi(numel(others)));
                end

                % 4) Compute the mixed payoff of focal node i using the original CDP payoff functions
                pair_avg_i   = pair_payoff_avg_CDP(i, old, Nhood, T, S, alpha, beta);
                triple_avg_i = triple_payoff_avg_CDP(i, old, nodeEdges, edgeNodes, T, S, G, W, alpha, beta);
                ben_i = (1 - delta) * pair_avg_i + delta * triple_avg_i;

                % 5) Compute the mixed payoff of role model j
                pair_avg_j   = pair_payoff_avg_CDP(j, old, Nhood, T, S, alpha, beta);
                triple_avg_j = triple_payoff_avg_CDP(j, old, nodeEdges, edgeNodes, T, S, G, W, alpha, beta);
                ben_j = (1 - delta) * pair_avg_j + delta * triple_avg_j;

                % 6) Fermi imitation: focal node i imitates role model j
                p_imitate = 1 / (1 + exp(-omega * (ben_j - ben_i)));
                if rand < p_imitate
                    state(i) = old(j);
                end

                % 7) Early stopping based on the stability of the C+P share
                if mod(t, checkEvery) == 0
                    cp_share = (sum(state==2) + sum(state==7)) / n;

                    if ptr < bufLen
                        ptr = ptr + 1;
                        cbuf(ptr) = cp_share;
                    else
                        cbuf(1:end-1) = cbuf(2:end);
                        cbuf(end)     = cp_share;
                    end

                    cc = cbuf(~isnan(cbuf));
                    if numel(cc) >= bufLen && std(cc) < tol_std
                        break;
                    end
                end
            end

            parCend(r) = sum(state==2) / n;
            parPend(r) = sum(state==7) / n;
            toc
        end

        % Store all repeated simulation results
        C_end(:, compIdx, dIdx)  = parCend;
        P_end(:, compIdx, dIdx)  = parPend;
        CP_end(:, compIdx, dIdx) = parCend + parPend;

        % Store mean values
        meanC(compIdx, dIdx)  = mean(parCend);
        meanP(compIdx, dIdx)  = mean(parPend);
        meanCP(compIdx, dIdx) = mean(parCend + parPend);

        fprintf('comp %d/%d finished at delta = %.2f | mean(C+P)=%.4f\n', ...
            compIdx, nComps, delta, meanCP(compIdx, dIdx));
    end
end

%% ===================== Save raw data =====================
save('CDP_equalProbNeighbor_MC_raw.mat', ...
    'C_end', 'P_end', 'CP_end', ...
    'meanC', 'meanP', 'meanCP', ...
    'comps', 'deltas');

%% ===================== Simple diagnostic plot based on mean values =====================
[~, midIdx] = min(abs(deltas - 0.5));
mc = meanCP(:, midIdx);

figure('Color','w'); hold on;
scatter3(comps(:,1), comps(:,2), comps(:,3), 50, mc, 'filled');
xlabel('Initial C'); ylabel('Initial D'); zlabel('Initial P');
title(sprintf('Mean C+P at \\delta = %.2f', deltas(midIdx)));
colorbar; box on; grid on; view(135,30);

%% ===================== Local functions =====================

function val = pair_payoff_avg_CDP(i, old, Nhood, T, S, alpha, beta)
    nei_i = Nhood{i};
    deg_i = numel(nei_i);
    si = old(i);

    if deg_i == 0
        val = 0;
        return;
    end

    pay_sum = 0;
    for u = 1:deg_i
        sj = old(nei_i(u));
        pay_sum = pay_sum + pay2_3str(si, sj, T, S, alpha, beta);
    end
    val = pay_sum / deg_i;
end

function val = triple_payoff_avg_CDP(i, old, nodeEdges, edgeNodes, T, S, G, W, alpha, beta)
    Hi = nodeEdges{i};
    m_i = numel(Hi);
    si = old(i);

    if m_i == 0
        val = 0;
        return;
    end

    acc = 0;
    for hh = 1:m_i
        v = edgeNodes{Hi(hh)};
        if numel(v) ~= 3
            continue;
        end
        s = [old(v(1)), old(v(2)), old(v(3))];
        acc = acc + pay3_3str(si, s, T, S, G, W, alpha, beta);
    end
    val = acc / m_i;
end

function val = pay2_3str(focal, other, T, S, alpha, beta)
    if focal == 0        % D
        if other == 0
            val = 0;
        elseif other == 2
            val = T;
        elseif other == 7
            val = T - beta;
        else
            val = 0;
        end
    elseif focal == 2    % C
        if other == 0
            val = S;
        else
            val = 1;
        end
    else                 % P
        if other == 0
            val = S - alpha;
        else
            val = 1;
        end
    end
end

function val = pay3_3str(focal, s, T, S, G, W, alpha, beta)
    nD = sum(s==0);
    nC = sum(s==2);
    nP = sum(s==7);

    switch focal
        case 0  % D
            if nD==3
                val = 0;
            elseif nC==1 && nD==2 && nP==0
                val = W;
            elseif nC==2 && nD==1 && nP==0
                val = T;
            elseif nP==1 && nD==2 && nC==0
                val = W - beta;
            elseif nP==1 && nC==1 && nD==1
                val = T - beta;
            elseif nP==2 && nD==1 && nC==0
                val = T - 2*beta;
            else
                val = 0;
            end

        case 2  % C
            if nC==1 && nD==2 && nP==0
                val = S;
            elseif nC>=1 && nD==1 && nP==0
                val = G;
            elseif nC==3 && nD==0 && nP==0
                val = 1;
            elseif nC>=1 && nD==1 && nP==1
                val = G;
            elseif nC>=1 && nD==0 && nP>=1
                val = 1;
            else
                val = 1;
            end

        case 7  % P
            if nP==1 && nD==2 && nC==0
                val = S - 2*alpha;
            elseif nP==1 && nD==1 && nC==1
                val = G - alpha;
            elseif nP==1 && nD==0 && nC==2
                val = 1;
            elseif nP==3 && nD==0 && nC==0
                val = 1;
            elseif nP==2 && nD==1 && nC==0
                val = G - alpha;
            elseif nP==2 && nD==0 && nC==1
                val = 1;
            else
                val = 1;
            end
    end
end