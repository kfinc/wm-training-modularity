%----------------------------------------------------
% Working memory training: Static normalized modularity calculation

% Karolina Finc | Centre for Modern Interdisciplinary Technologies, Nicolaus Copernicus University in Toruń, Poland
% Last edited: 09-01-2019
%----------------------------------------------------

% Mere the input is 5D matrix where:
% 1st dim --> subject
% 2nd dim --> session
% 3rd dim --> condition
% 4th and 5th dim --> correlation matrix

%% Loading data

M = load('WM_all_static_power_correlation_matrices.mat');

M = M.all_conditions_new_power;
M = squeeze(M);

%% Params
n_sub = size(M,1);
n_ses = size(M,2);
n_con = size(M,3);
n_rep = 100;
n_null = 100;
n_rew = 1;

q = zeros(n_sub, n_ses, n_con);
q_null_mean = zeros(n_sub, n_ses, n_con);

%% Static modularity calculation
for sub = 1: n_sub
    sub
    for ses = 1: n_ses
        for con = 1: n_con
            Qb = 0;
            for rep = 1: n_rep
                A = squeeze(M(sub,ses,con,:,:)).* squeeze(M(sub, ses, con, :, :) > 0);
                [~, Qt] = community_louvain(A, 1, []);
                 if Qt > Qb 
                     Qb = Qt;
                 end
                 q(sub, ses, con) = Qb;
            end  
    end
    end

end
save('mean_power_static_modularity.mat', 'q');

%% Randomized modularity calculation
for sub = 1: n_sub
    q_null = zeros(1, n_null);

    sub
    for ses = 1: n_ses
        for con = 1: n_con
            A = squeeze(M(sub,ses,con,:,:)).* squeeze(M(sub, ses, con, :, :) > 0);
            if sum(isnan(A(:))) > 0
                q_null = 0;
            else    
                for null = 1 : n_null
                    B = randmio_und(A, n_rew);
                    Qb = 0;        
                    for rep = 1: n_rep 
                        [~, Qt] = community_louvain(B, 1, []);
                        if Qt > Qb 
                            Qb = Qt;
                        end
                    end
                    q_null(1, null) = Qb;
                end
            end
           q_null_mean(sub, ses, con) = mean(q_null);

        end
    end
end

%% Normalized modularity calculation
q_norm = q_mean ./ q_null_mean;
save('mean_power_normalized_modularity.mat', 'q_norm');