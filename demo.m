close all;
% clc;
clear;
Nv = 30;
Ne = 50;
Ls = 1000;
max_iter = 100;
%% Generating undirected graph, Adjacency & Laplacian
A = rand_ugraph(Nv,Ne,0.1,0.4);
L = diag(sum(A)) - A;

%% Generating test signals
%[V, Lbd] = eig(L);
%sp = rand(2,Ls);
%sp = [sp; zeros(Nv - 2, Ls)];
%S = V*sp;
%mu = ones(Nv,1);
%sigma = pinv(Lbd);
%R = chol(sigma);
%R = sigma;
%sp = repmat(mu, 1, Ls) + (randn(Ls, Nv)*R)'; 
%S = V*sp;
% Adding random normal distributed noise
%S = S + randn(Nv, Ls)./10;
stmls = randn(Nv, Ls);
S = (A - eye(Nv))\stmls;
%% Calling GL_SigRep to estimate Laplacian
disp('SigRep with quadprog')
tic
[Sr, Le] = GL_SigRep(S, 0.1, 20, max_iter);
toc
disp('SigRep with CVX');
tic
[Sr_CVX, Le_CVX] = GL_SigRep_CVX(S, 0.1, 20, max_iter);
toc
%% Ploting
close all;

figure;
subplot(2, 3, 1)
G = graph(A);
pG = plot(G);
pG.LineWidth = 5*G.Edges.Weight/max(G.Edges.Weight);
title('Raw Graph');
subplot(2, 3, 4)
imagesc(L);
title('Raw Laplacian');
subplot(2, 3, 2)
Ae = -(Le-diag(diag(Le)));
Ae(abs(Ae) < 1e-5) = 0;
Ge = graph(Ae);
pGe = plot(Ge);
pGe.XData = pG.XData;
pGe.YData = pG.YData;
pGe.LineWidth = 5*Ge.Edges.Weight/max(Ge.Edges.Weight);
title('Estimated Graph with quadprog');
subplot(2, 3, 5)
imagesc(Le);
title('Estimated Laplacian with quadprog');


subplot(2, 3, 3)
Ae_CVX = -(Le_CVX-diag(diag(Le_CVX)));
Ae_CVX(abs(Ae_CVX) < 1e-5) = 0;
Ge_CVX = graph(Ae_CVX);
pGe_CVX = plot(Ge_CVX);
pGe_CVX.XData = pG.XData;
pGe_CVX.YData = pG.YData;
pGe_CVX.LineWidth = 5*Ge_CVX.Edges.Weight/max(Ge_CVX.Edges.Weight);
title('Estimated Graph with CVX');
subplot(2, 3, 6)
imagesc(Le_CVX);
title('Estimated Laplacian with CVX');