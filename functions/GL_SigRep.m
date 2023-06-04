function [Y, L] = GL_SigRep(X, alpha, beta, max_iter)
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
vYY = mat2vec(Y*Y');


%% Estimation Iterations

cvx_begin
    variable l((1+n)*n/2)
    minimise (alpha*vYY'*Md*l + beta*l'*Md'*Md*l)
    subject to
    A*l <= a
    B*l == b
cvx_end

L = vec2mat(Md*l, n);
%for i = 1:max_iter
   

%end
