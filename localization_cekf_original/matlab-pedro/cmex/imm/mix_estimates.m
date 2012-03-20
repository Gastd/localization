%Mixes the state vectors (old_state_vectors) and covariance matrices 
%(old_cov_matrices) from the previous iteration and outputs the mixed 
%state vectors (mixed_state_vectors) and covarince matrices 
%(mixed_cov_matrices) that will be used in the next filtering step. The 
%mixing uses the mode probability vectors from the previous iteration 
%(old_mode_prob_vector) and the predicted vectors (pred_mode_prob_vector) 
%according to the TPM (Markov_matrix)
function [mixed_state_vectors,mixed_cov_matrices] = mix_estimates(old_state_vectors,old_cov_matrices,old_mode_prob_vector,pred_mode_prob_vector,Markov_matrix)
    
    %Retrieves the number of vectors and their dimensions
    dimVectors = size(old_state_vectors);
    %Retrieves the number of matrices and their dimensions
    dimMatrices = size(old_cov_matrices);
    
    %Memory allocation
    mixed_state_vectors = zeros(dimVectors);
    mixed_cov_matrices = zeros(dimMatrices);
    
    %dimVectors(3) = number of state vectors
    for i=1:dimVectors(3)
        
        %Mixes the state vectors
        for j=1:dimVectors(3)
            mixed_state_vectors(:,:,i) = mixed_state_vectors(:,:,i) + Markov_matrix(i,j)*old_mode_prob_vector(j,1)*old_state_vectors(:,:,j);
        end
        
        mixed_state_vectors(:,:,i) = mixed_state_vectors(:,:,i)/pred_mode_prob_vector(i,1);
        
    end
    
    %dimVectors(3) = number of matrices
    for i=1:dimMatrices(3)
        
        %Mixes the covariance matrices
        for j=1:dimMatrices(3)
            
            %Mixes the covariance matrices
            mixed_cov_matrices(:,:,i) = mixed_cov_matrices(:,:,i) + ...
                Markov_matrix(i,j)*old_mode_prob_vector(j,1)*( old_cov_matrices(:,:,j) + ...
                (old_state_vectors(:,:,j)-mixed_state_vectors(:,:,i))*(old_state_vectors(:,:,j)-mixed_state_vectors(:,:,i))' );
        end
        
        mixed_cov_matrices(:,:,i) = mixed_cov_matrices(:,:,i)/pred_mode_prob_vector(i,1);
        
    end