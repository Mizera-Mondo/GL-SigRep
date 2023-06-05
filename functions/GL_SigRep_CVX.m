function [Y, L] = GL_SigRep_CVX(X, alpha, beta, max_iter)
%% Using GL_SigRep algorithm to estimate underlying 
% graph Laplacian and denoised real graph signal matrix
% So: Signal Observed, N-by-K matrix
% alpha: Regularization parametre of signal stationarity, higher for more stationary
%           denoised signal under estimated Laplacian 
% beta: Regularization parametre of sparsity of Laplacian
% max_iter: maximum iteration

% Sr: Signal Recovered(Denoised), N-by-K matrix 
% L: Laplacian estimated, N-by-N square matrix

%% Initialization
[n,~] = size(X);
Y = X;
Z = zeros(n);
I = eye(n);
Md = dupmat(n);

% Constraints
% Inequities, A*lh <= 0
A = diag(mat2vec(ones(n,n) - I))*dupmat(n);
[m,~] = size(A);
a = zeros(m,1);
% Equalties, B*lh - b = 0 
B = kron(ones(1,n),I); % L*1 = 0
B = [B; (mat2vec(I))']; % tr(L) = n
B = B*Md;
b = zeros(n+1,1);
b(n+1) = n;



%% Estimation Iterations
L = randn(n);

for i = 1:max_iter
    L_old = L;
    Y_old = Y;
    vYY = mat2vec(Y*Y');
    % Sub-problem of L with CVX
    cvx_begin quiet
        variable l((1+n)*n/2)
        minimise (alpha*vYY'*Md*l + beta*l'*Md'*Md*l)
        subject to
        A*l <= a
        B*l == b
    cvx_end
    L = vec2mat(Md*l, n);
    
    % Sub-problem of Y
    Y = (I + alpha.*L)\X;

    % Checking terminating condition
    err1 = norm(L - L_old, 'fro');
    err2 = norm(Y - Y_old, 'fro');
    if err1 < 1e-6 && err2 < 1e-6
        break;
    end
end
end
