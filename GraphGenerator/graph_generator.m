
function [Theta,Sigma]=graph_generator(p,seed,eigg, output_folder, varargin)
% Generate inverse covariance and covariance matrix under certain graph structures
% Inputs:
%%% p: total number of features.
%%% seed: random starting seed.
%%% eigg: low-rankness of generated inverse covariance matrix.
%%%% 1/eigg is the leading eigenvalue of the covariance matrix;
%%%% All other eigenvalues will be <= 1. 
%%% output_folder: directory to save results.
%%% varagin: type of graph to produce for the inverse covariance matrix.
% Outputs:
%%% Theta: inverse covariance matrix. 
%%% Sigma: covariance matrix. (Will not be sparse.)

% Initialization
Theta=zeros(p,p);
rng(seed);
type=varargin;
mkdir(output_folder);
% Simple chain graph
if type=="chain"
    disp("here")
    for i=1:p
        if i<p
            Theta(i,i+1)=(rand(1)*2-1);
            Theta(i,i+1)=Theta(i,i+1)+sign(Theta(i,i+1));
        end
        Theta(i,i)=rand(1)+1;
        if i>1
            Theta(i,i-1)=Theta(i-1,i);
        end
    end
% Block diagonal graph.    
elseif type=="block"
    blocksize= floor(p / 3);
    K=floor(p/blocksize);
    prob=1; 
    a=rand([p,p]);
    edges= triu(a<prob);
    edges = edges + edges';
    for i=1:K
        if i<K
            set=((i-1)*blocksize+1):(i*blocksize);
        else
            set=((i-1)*blocksize+1):p;
        end
        for j=1:length(set)
            Theta(set(j),set(j))=rand(1)+1;
            Theta(set(j),set((j+1):length(set)))=rand(1)*2-1;
            Theta(set(j),set((j+1):length(set)))=Theta(set(j),set((j+1):length(set)))...
                +sign(Theta(set(j),set((j+1):length(set))));
            Theta(set(j), :)=Theta(set(j), :).*edges(set(j), :);
            Theta(set((j+1):length(set)),set(j))=Theta(set(j),set((j+1):length(set)))';
        end 
    end
% Erdos-Renyi graph.     
elseif type=="ER"
    prob=0.03;
    a=rand([p,p]);
    edges=tril(a<prob);
    Theta= Theta + rand([p,p])* 0.45 + 0.15;
    Theta=Theta+sign(Theta);
    Theta=Theta.*edges;
    Theta=Theta+Theta';
    for j=1:p
        Theta(j,j)=rand(1);
        Theta(j,j)=Theta(j,j)+sign(Theta(j,j));
    end
% Simple star graph
elseif type == "star"
    for j=1:p
        Theta(j,j)=rand(1);
        Theta(j,j)=Theta(j,j)+sign(Theta(j,j));
        if j>1
            Theta(1,j)=rand(1)*2-1;
            Theta(1,j)=Theta(1,j)+sign(Theta(1,j));
            Theta(j,1)=Theta(1,j);
        end
    end
% Watts-Strogatz small-world graph. 
elseif type == "WS"
    k= round(p * 1.5) ; prob=0.1;
    %2*k is the nsaumber of neighbors
    %prob is the probability that a neighbor is substituted by a randomly
    %picked one
    g=WattsStrogatz(p,k,prob);
    A=adjacency(g);
    Theta=diag(rand(1,p))/2;
    rr = (rand(p,nnz(A)/2)*2-1);
    Theta(tril(A~=0)) = rr(tril(A~=0));
    Theta=Theta+Theta';
    Theta=Theta+sign(Theta);
else 
    disp('Non-valid graph type input; defaulting to ER graph.')
    prob=0.03;
    a=rand([p,p]);
    edges=tril(a<prob);
    Theta= Theta + rand([p,p])* 0.45 + 0.15;Theta=Theta+sign(Theta);
    Theta=Theta.*edges;
    Theta=Theta+Theta';
    for j=1:p
        Theta(j,j)=rand(1);
        Theta(j,j)=Theta(j,j)+sign(Theta(j,j));
    end
end
% Implement low-rank covariance matrix.
Theta=Theta+(eigg-min(eig(Theta)))*eye(p);
Sigma = inv(Theta);
% Save results.
writematrix(Theta, strcat(output_folder, 'TruePrecision.csv'));
writematrix(Sigma, strcat(output_folder, 'FullCov.csv'));
end
