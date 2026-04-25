function out = NOTEARS(X)
% NOTEARS
%
% INPUT
%   X : [T x d] data matrix (rows = samples, columns = variables)
%
% OUTPUT (struct out)
%   out.connectivity : [d x d] weighted adjacency matrix, out.connectivity(i,j)=i->j
%   out.flows        : [E x 1] edge weights aligned with out.edges (same ordering)
%   out.edges        : [E x 2] directed edge list [i j] for nonzero edges in connectivity
%   out.h_final      : final acyclicity value h(W)
%
% Reference:
%  Zheng, X., Bryon A., Pradeep K. R., and Eric P. X. 2018 Dags with no tears: 
%  continuous optimization for structure learning. Advances in 
%  Neural Information Processing Systems 31.
%
% (C) 2026 Moo K. Chung

[T,d] = size(X);
X = zscore(X);

% ----- hyperparameters -----
lambda1 = 0.01;     % L1 sparsity
rho     = 1.0;      % augmented Lagrangian penalty
alpha   = 0.0;      % Lagrange multiplier
hTol    = 1e-8;     % acyclicity tolerance
maxIter = 50;
wThresh = 1e-3;

nInner  = 50;

% ----- precompute Gram -----
S = (X' * X) / T;
L0 = norm(S,2);
if L0 == 0, L0 = 1; end

% ----- initialize -----
W = zeros(d,d);
h_prev = inf;

for iter = 1:maxIter

    for s = 1:nInner

        % loss gradient for 0.5*||X - XW||_F^2 / T  (up to constants)
        G_loss = S*W - S;

        % acyclicity h(W) and gradient
        E  = expm(W.*W);
        h  = trace(E) - d;
        G_h = (E') .* (2*W);

        % augmented Lagrangian gradient
        G = G_loss + (rho*h + alpha)*G_h;

        step = 1 / (L0 + rho*max(1, norm(G_h,'fro')) + 1e-6);

        W = soft_threshold(W - step*G, step*lambda1);
        W(1:d+1:end) = 0;
    end

    E = expm(W.*W);
    h = trace(E) - d;

    if abs(h) <= hTol
        break
    end

    alpha = alpha + rho*h;

    if abs(h) > 0.25*abs(h_prev)
        rho = min(rho*2.0, 1e6);
    end
    h_prev = h;
end

W(abs(W) < wThresh) = 0;
W(1:d+1:end) = 0;

% ----- baseline-style outputs -----
out = struct();
out.connectivity = W;
out.h_final      = h;

[ei, ej, ev] = find(W);
out.edges = [ei ej];
out.flows = ev;

end

% --------------------------------------------------------------
function Y = soft_threshold(X, tau)
Y = sign(X) .* max(abs(X) - tau, 0);
end
