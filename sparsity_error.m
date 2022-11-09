function mean_errors=sparsity_error(alphas, n1, n2, q, samplesize)
% Computes mean error for a stochastic block model 
%           with 2 communities and connectivity matrix
%           alpha*B, where B=[1 q;  
%                             q 1]
%           for values of alpha given as function input. For each
%           alpha a sample of sample size specified as input
%           is drawn from the model.
%
%     Input parameters
%       alphas - vector of scaling parameters alpha
%       n1 - size of the first community
%       n2 - size of the second community
%       q - edge probability for vertices of different communities in 
%           the unscaled connectivity matrix
%       samplesize - number of adjacency matrices that are sampled for
%           each scaling parameter
%
%     Output value
%       mean_errors - an array containing the mean fraction of
%           misclassified vertices of the sampled adjacency matrices 
%           for each scaling parameter in alphas
%
    %% Set parameters
    % randon number generator for reproducable results
    rng('default');
    % compute total number of vertices
    n = n1 + n2;

    % vector for community assignment
    Assignment = zeros(n,1);
    % first n1 vertices are assigned to community 1
    Assignment(1:n1) = 1;
    % vertices n1+1 through n are assigned to community 2
    Assignment(n1+1:n) = 2;

    % connectivity matrix B containing the inter-community edge
    % probabilities
    B= [1 q;
        q 1];

    % vector collecting the relative clustering errors
    % one row consists of the sample outcomes for one alpha
    relative_errors_total = zeros(length(alphas),samplesize);


    % counter for indexing
    counter = 0;
    
    % go through vector of scaling parameters
    for alpha=alphas
        
        counter = counter + 1;

        % scale B by the scaling parameter alpha
        B_scaled = alpha*B;

        % vector collecting the relative errors for each sampled
        % adjacency matrix
        rel_errors_alpha = zeros(1, samplesize);

        % sampling of adjacency matrix and error computation
        for sample_index=1:samplesize
            %% sample adjacency matrix
            % initialize an nxn matrix
            A=zeros(n,n);

            % sample entries on and above the diagonal as Bernoulli(p)-
            % distributed random variables with p given by the 
            % connectivity matrix B
            for i=1:n
                for j=i:n
                    A(i,j) = random('Binomial',1, ...
                        B_scaled(Assignment(i),Assignment(j)));

                    % set entries below the diagonal such that a 
                    % symmetric matrix is obtained
                    if i ~= j
                        A(j,i) = A(i,j);
                    end
                end
            end

            %% Run spectral clustering algorithm
            % compute 2 largest absolute eigenvalues and corresponding
            % eigenvectors
            [V,D] = eigs(A,2);

            % classify vertices based on k-means applied to the matrix 
            % with columns consisting of the 2 eigenvectors
            classification = kmeans(V,2);

            %% error is smallest errors over possible labelings
            % labeling 1
            error_vector_1 = (Assignment == classification);

            % labeling 2
            error_vector_2 = (Assignment ~= classification);

            % absolute error is minimum of absolute errors
            % over 2 possible labelings
            abs_error = min(sum(error_vector_1), sum(error_vector_2));

            % relative error as fraction of absolute misclassified 
            % vertices
            rel_error = abs_error/n;

            % add error for sampled adjacency matrix to vector
            % for current alpha
            rel_errors_alpha(sample_index) = rel_error;
        end
        % add relative errors obtained for current alpha to matrix 
        % collecting the errors for all alphas
        relative_errors_total(counter,:) = rel_errors_alpha;

    end

    % compute the mean fraction of misclassified vertices 
    % over all sampled adjacency matrices
    mean_errors = mean(relative_errors_total,2);
end