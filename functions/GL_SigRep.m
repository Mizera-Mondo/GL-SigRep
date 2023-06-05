function [Sr, L] = GL_SigRep(So, alpha, beta, max_iter)
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
[n,~] = size(So);
Sr = So;
L = zeros(n,n);
lh = sym2vech(L);
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

% Target function parametres
p1 = alpha*(mat2vec(Sr*Sr'))'*Md;
p2 = beta*(Md'*Md);

%% Estimation Iterations
for i = 1:max_iter

    % Calculating L, using Interior Point Method(IPM)
    % Or using quadprog()
    lh = quadprog(2*p2, p1, A, a, B, b, [], [] ,[], optimoptions('quadprog', 'Display', 'none'));
    Lnew = vec2squ(Md*lh);
    % Calculating Y
    Sr_new = (I + alpha.*Lnew)\So;
    %Updating p1 parametre
    p1 = alpha*(mat2vec(Sr_new*Sr_new'))'*Md;

    % Checking terminating condition
    err1 = norm(Lnew - L, 'fro');
    err2 = norm(Sr_new - Sr, 'fro');
    Sr = Sr_new;
    L = Lnew;
    if err1 < 1e-6 && err2 < 1e-6
        break;
    end
end
end

