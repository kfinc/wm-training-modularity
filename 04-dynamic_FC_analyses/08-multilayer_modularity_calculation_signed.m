clc; clear;

%% Loading data
large = load('LB_dualnback_power_dynamic_correlation_matrices.mat');
M = squeeze(large.correlation_matrices_dyn_wei);

n_sub = size(M,1);
n_ses = size(M,2);
T = size(M,3);
N = size(M,4);

%% Parameters
gplus = 1;
gminus = 1;
omega = 1;
n_rep = 100;

%% Empty matrices to store results
modularity_mean = zeros(n_sub, n_rep, n_ses); % Mean modularity
modules = zeros(n_sub, n_ses, n_rep, N, T); % Module assigment labels    

for sub = 1 : 1:n_sub
    for ses = 1 :1: n_ses
        
    %--- define objects ---------------------------------------------------
        A = cell(1, T);
        B=spalloc(N*T,N*T,N*N*T+2*N*T);
        twom=0;
    
        % Signed version by Xiaosong He, 5/23/2019
        for s=1:T
            A{s} = squeeze(M(sub, ses, s, :, :));
            Aplus=A{s}; Aplus(A{s}<0)=0;
            Aminus=-A{s}; Aminus(A{s}>0)=0;
            kplus=sum(Aplus)';
            kminus=sum(Aminus)';
            mm=sum(kplus)+sum(kminus); 
            twom=twom+mm;
            indx=[1:N]+(s-1)*N;
        if sum(kminus) == 0 % necessary for case no negative edges in matrix
            Gm=0;
        else
            Gm=gminus*kminus*kminus'/sum(kminus);
        end
            B(indx,indx)=A{s}-(gplus*kplus*kplus'/sum(kplus)-Gm);
        end
        

        B = B + omega/2*spdiags(ones(N*T,2),[-N, N], N*T, N*T);
        B = B + omega*spdiags(ones(N*T,2),[-2*N, 2*N], N*T, N*T);

        twom=twom+(N*T*(T-1)*omega);
        
        for rep = 1 : n_rep
            clc;
            fprintf('Subject = %i\n',sub);
            [S,Q] = genlouvain(B);
            Q = Q / twom;
            S = reshape(S, N, T);
            
            modularity_mean(sub, rep, ses) = Q;
            modules(sub, ses, rep, :, :) = S;
        end
        
   end
end
