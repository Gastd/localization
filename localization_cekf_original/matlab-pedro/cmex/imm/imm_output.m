%Gives the IMM Algorithm's output based on the state estimates given by the
%multiple estimators
function [xk_est_IMM,Pk_est_IMM] = imm_output(corr_state_vectors,corr_cov_matrices,updated_mode_prob_vector)

    %Retrieves the number of vectors and their dimensions
    dimVectors = size(corr_state_vectors);
    %Retrieves the number of matrices and their dimensions
    dimMatrices = size(corr_cov_matrices);
    
    %Memory allocation
    xk_est_IMM = zeros(dimVectors(1),dimVectors(2));
    Pk_est_IMM = zeros(dimMatrices(1),dimMatrices(2));
    
     %dimVectors(3) = number of state vectors
    for i=1:dimVectors(3)
        xk_est_IMM = xk_est_IMM + updated_mode_prob_vector(i,1)*corr_state_vectors(:,:,i);        
    end
    
    for i=1:dimVectors(3)
            Pk_est_IMM = Pk_est_IMM + updated_mode_prob_vector(i,1)*(corr_cov_matrices(:,:,i)+(corr_state_vectors(:,:,i)-xk_est_IMM)*(corr_state_vectors(:,:,i)-xk_est_IMM)');
    end