%% ===================== Mean C vs. e and delta (equal-probability 2/3 neighbor selection) =====================
% Inputs in workspace:
%   B : (n x n) 0/1, undirected, no self loops
%   M : (n x H) each column is a 3-uniform hyperedge (exactly 3 ones)
%
% Description:
% 1) At each step, a focal node i is randomly selected.
% 2) If node i has both pairwise and three-body neighbors, the update channel is selected with equal probability.
% 3) If only one type of neighbor exists, that channel is used.
% 4) Payoffs are still computed using the original pair_avg and triple_avg terms, then weighted by delta.
% 5) The focal node i updates its strategy according to the Fermi rule.

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

rng(42);   % reproducible

%% -------------------- Basic checks --------------------
B = spones(B);
n = size(B,1);

% Force undirected structure and remove self-loops
B = triu(B,1) + triu(B,1)';
assert(isequal(B,B.'), 'B must be symmetric, corresponding to an undirected graph.');
assert(all(diag(B)==0), 'B must not contain self-loops.');

% Ensure 3-uniform hyperedges
M = spones(M);
if size(M,1) ~= n
    error('The number of rows in M must be equal to the size of B.');
end

colsum = full(sum(M,1));
keep = (colsum == 3);
if any(~keep)
    warning('Removed %d non-3-uniform hyperedge columns.', nnz(~keep));
    M = M(:, keep);
end
H = size(M,2);

fprintf('n = %d\n', n);
fprintf('2-body edges = %d\n', nnz(triu(B,1)));
fprintf('3-body hyperedges = %d\n', H);

%% -------------------- Precompute neighborhoods / hyperedges --------------------
Nhood = cell(n,1);        % Pairwise neighbors
for i = 1:n
    Nhood{i} = find(B(i,:));
end

edgeNodes = cell(H,1);    % Nodes in each 3-hyperedge
for h = 1:H
    edgeNodes{h} = find(M(:,h))';   % 1 x 3
end

nodeEdges = cell(n,1);    % Incident 3-hyperedges of each node
for i = 1:n
    nodeEdges{i} = find(M(i,:));
end

deg2 = cellfun(@numel, Nhood);
deg3 = cellfun(@numel, nodeEdges);

fprintf('nodes with 2-body neighbors: %d / %d\n', nnz(deg2>0), n);
fprintf('nodes with 3-body neighbors: %d / %d\n', nnz(deg3>0), n);

%% -------------------- Payoff / dynamics parameters --------------------
omega = 0.1;      % Fermi intensity

% Payoff parameters
S = -0.1;
T = 1.8;
G = 0.8;
W = 0.1;

% Delta only represents the weight of the three-body payoff contribution
deltas = 0.0:0.1:1.0;
eGrid  = 0.05:0.05:0.95;

rep  = 100;
tmax = 500000;

% Early stopping
checkEvery = 500;
bufLen     = 10;
tol_std    = 1e-4;

%% -------------------- Result containers --------------------
E = numel(eGrid);
D = numel(deltas);

meanC = zeros(E, D);
C_end = zeros(rep, E, D);

%% ===================== Main loops =====================
for dIdx = 1:D
    delta = deltas(dIdx);

    fprintf('\n===== delta = %.2f =====\n', delta);

    for eIdx = 1:E
        eC = eGrid(eIdx);
        nC = round(n * eC);
        nC = max(0, min(n, nC));

        parCend = zeros(rep,1);

        parfor r = 1:rep
            % ---------- Initialization: D = 0, C = 2 ----------
            state = zeros(n,1,'int8');    % All D
            if nC > 0
                idxC = randperm(n, nC);
                state(idxC) = 2;          % C
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
                    v = edgeNodes{h};   % 1 x 3
                    others = v(v ~= i);
                    j = others(randi(numel(others)));
                end

                % 4) Compute the mixed payoff of focal node i using the original formula
                pair_avg_i   = pair_payoff_avg(i, old, Nhood, T, S);
                triple_avg_i = triple_payoff_avg(i, old, nodeEdges, edgeNodes, T, S, G, W);
                ben_i = (1 - delta) * pair_avg_i + delta * triple_avg_i;

                % 5) Compute the mixed payoff of role model j using the original formula
                pair_avg_j   = pair_payoff_avg(j, old, Nhood, T, S);
                triple_avg_j = triple_payoff_avg(j, old, nodeEdges, edgeNodes, T, S, G, W);
                ben_j = (1 - delta) * pair_avg_j + delta * triple_avg_j;

                % 6) Fermi adoption: focal node i imitates role model j
                p = 1 / (1 + exp(-omega * (ben_j - ben_i)));
                if rand < p
                    state(i) = old(j);
                end

                % 7) Early stopping based on the stability of the C share
                if mod(t, checkEvery) == 0
                    cshare = sum(state==2) / n;

                    if ptr < bufLen
                        ptr = ptr + 1;
                        cbuf(ptr) = cshare;
                    else
                        cbuf(1:end-1) = cbuf(2:end);
                        cbuf(end) = cshare;
                    end

                    cc = cbuf(~isnan(cbuf));
                    if numel(cc) >= bufLen && std(cc) < tol_std
                        break;
                    end
                end
            end

            parCend(r) = sum(state==2) / n;
        end

        C_end(:, eIdx, dIdx) = parCend;
        meanC(eIdx, dIdx) = mean(parCend);

        fprintf('e = %.2f, mean terminal C = %.4f\n', eC, meanC(eIdx, dIdx));
    end
end

%% ===== Save raw data =====
save('CD_equalProbNeighbor_MC_raw.mat', 'C_end', 'meanC', 'eGrid', 'deltas');

%% ===== Visualization =====
figure('Color','w');
imagesc(deltas, eGrid, meanC);
axis xy
colormap(cividis);
colorbar;
xlabel('\delta');
ylabel('Initial C share, e');
title('Mean terminal proportion of C (equal-probability 2/3 neighbor selection)');

figure('Color','w');
hold on
co = [31 119 180;
      255 127 14;
      148 103 189;
      44 160 44;
      214 39 40;
      23 190 207] / 255;

for dIdx = 1:D
    plot(eGrid, meanC(:,dIdx), '-', 'LineWidth', 2, ...
        'Color', co(1+mod(dIdx-1,size(co,1)),:));
end
xlabel('Initial C share, e');
ylabel('Mean terminal C share');
box on; grid on
legend(arrayfun(@(x) sprintf('\\delta = %.1f', x), deltas, 'UniformOutput', false), ...
       'Location','northwest','Box','off');

%% ===================== Local functions =====================

function val = pair_payoff_avg(i, old, Nhood, T, S)
    nei_i = Nhood{i};
    deg_i = numel(nei_i);

    if deg_i > 0
        c_cnt = sum(old(nei_i)==2);
    else
        c_cnt = 0;
    end

    if old(i)==0      % D
        val = (deg_i>0) * (T * c_cnt / max(deg_i,1));
    else              % C
        val = (deg_i>0) * ((1 * c_cnt + S * (deg_i - c_cnt)) / max(deg_i,1));
    end
end

function val = triple_payoff_avg(i, old, nodeEdges, edgeNodes, T, S, G, W)
    Hi = nodeEdges{i};
    m_i = numel(Hi);

    if m_i == 0
        val = 0;
        return;
    end

    acc = 0;
    for hh = 1:m_i
        v = edgeNodes{Hi(hh)};   % 1 x 3
        ssum = double(old(v(1))) + double(old(v(2))) + double(old(v(3)));
        acc = acc + pay3_2str(old(i), ssum, T, S, G, W);
    end
    val = acc / m_i;
end

function val = pay3_2str(focal, ssum, T, S, G, W)
% States: D = 0, C = 2
% ssum belongs to {0, 2, 4, 6}
    if focal==0      % D
        switch ssum
            case 0
                val = 0;    % D-D-D
            case 2
                val = W;    % C-D-D
            case 4
                val = T;    % C-C-D
            case 6
                val = T;    % C-C-C
            otherwise
                val = 0;
        end
    else             % C
        switch ssum
            case 0
                val = S;    % Rare case
            case 2
                val = S;    % C-D-D
            case 4
                val = G;    % C-C-D
            case 6
                val = 1;    % C-C-C
            otherwise
                val = 1;
        end
    end
end