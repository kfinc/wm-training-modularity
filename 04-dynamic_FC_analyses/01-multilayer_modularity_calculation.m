%----------------------------------------------------
% Multilayer community detection algoritm

% Karolina Finc | Centre for Modern Interdisciplinary Technologies, Nicolaus Copernicus University in Toruń, Poland
% Last edited: 12-09-2018
%----------------------------------------------------

% Mere the input is 5D matrix where:
% 1st dim --> subject
% 2nd dim --> session
% 3rd dim --> time windows
% 4th and 5th dim --> correlation matrix

clc; clear;

%% Loading data
large = load('LB_dualnback_power_dynamic_correlation_matrices.mat');
M = squeeze(large.correlation_matrices_dyn_wei);

n_sub = size(M,1);
n_ses = size(M,2);
n_win = size(M,3);
n_roi = size(M,4);


%% Parameters
gamma = 1;
omega = 1;
n_rep = 100;

%% Empty matrices to store data
modularity_mean = zeros(n_sub, n_rep, n_ses); % Mean modularity
modules = zeros(n_sub, n_ses, n_rep, n_roi, n_win); % Module assigment labels    

%%
for sub = 1 : n_sub
    for ses = 1 : n_ses
    %--- define objects ---------------------------------------------------
        A = cell(1, n_win);
        B = spalloc(n_roi * n_win, n_roi * n_win,(n_roi + n_win) * n_roi* n_win);
        twomu = 0;
    
        %--- null model -------------------------------------------------------
        for win = 1 : n_win
            %--- copy network with positive weights thresholding --------------
            A{win} = squeeze(M(sub, ses, win, :, :) .* (M(sub, ses, win, :, :) > 0));
            k = sum(A{win});                             % node degree
            twom = sum(k);                               % mean network degree
            twomu = twomu + twom;                        % increment
            indx = [1:n_roi] + (win-1)*n_roi;            % find indices
            B(indx,indx) = A{win} - gamma * [k'*k]/twom;  % fill B matrix
        end
        twomu = twomu + 2*omega* n_roi*(n_win-1);

        B = B + omega/2*spdiags(ones(n_roi*n_win,2),[-n_roi, n_roi], n_roi*n_win, n_roi*n_win);
        B = B + omega*spdiags(ones(n_roi*n_win,2),[-2*n_roi, 2*n_roi], n_roi*n_win, n_roi*n_win);
        %B = B + omega * spdiags(ones(n_roi * n_win,2),[-n_roi, n_roi], n_roi*n_win, n_roi*n_win);

%--- calculate multilayer modules -------------------------------------
        %Qb = 0;
        for rep = 1 : n_rep
            clc;
            fprintf('Subject = %i\n',sub);
            [S,Q] = genlouvain(B);
            Q = Q / twomu;
            S = reshape(S, n_roi, n_win);
            
            modularity_mean(sub, rep, ses) = Q;
            modules(sub, ses, rep, :, :) = S;
       end
    end
end
