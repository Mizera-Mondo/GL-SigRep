close all;
clc;
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
[Sr, Le] = GL_SigRep(S,0.1,10,max_iter);

%% Ploting
close all;
figure;
G = graph(A);
pG = plot(G);
pG.LineWidth = 5*G.Edges.Weight/max(G.Edges.Weight);
title('Raw Graph');
figure;
imagesc(L);
title('Raw Laplacian');
figure;
Ae = -(Le-diag(diag(Le)));
Ae(abs(Ae) < 1e-5) = 0;
Ge = graph(Ae);
pGe = plot(Ge);
pGe.XData = pG.XData;
pGe.YData = pG.YData;
pGe.LineWidth = 5*Ge.Edges.Weight/max(Ge.Edges.Weight);
title('Estimated Graph');
figure;
imagesc(Le);
title('Estimated Laplacian');