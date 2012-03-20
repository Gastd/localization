%   Transition Probability Matrix (TPM) estimator for Markov Chains
%   following an a priori Dirichlet distribution. This simulation does not
%   assume any mode knowledge. It uses just the measurements' likelihoods
%   with respect to each possible mode. This code is based on the work
%
%   "Online Bayesian Estimation of Transition Probabilities for Markovian
%   Jump Systems"
% 
%   Author: Pedro Henrique R.Q.A. Santana
%
%   Parameters:
%
%       yk -> measurement
%       xk_pred -> predicted state vector
%       Pk_pred -> covariance for the predicted state vector
%       C_matrices -> set of output matrices for each mode
%       R_matrices -> set of covariances matrices for each mode's measurement
%       TPM_ant -> previous estimative of the TPM
%       prob_vector_ant -> 
%       alpha_matrix_ant -> previous values for the Dirichlet Distributions' alpha parameters 
%       iteration_number -> iteration counter (starts at 0)
%       inital_alpha_sum -> initial sum of the alpha parameters for each column
%
%       TPM(i,j) = P{xk+1=i|xk=j} -> The columns add to one!
%


function [TPM_pos,alpha_matrix_pos] = likelihoodDirichletEstimator(likelihood_vect,TPM_ant,prob_vector_ant,alpha_matrix_ant,iteration_number,inital_alpha_sum)

%   Takes the matrix transpose to make it organizes in rows.
TPM_ant = TPM_ant';
alpha_matrix_ant = alpha_matrix_ant';

%   Probability vector's dimension
dimProbVect = size(prob_vector_ant);

alpha_matrix_pos = alpha_matrix_ant;    %   Memory allocation
g_matrix = TPM_ant;                     %   Memory allocation
TPM_pos = TPM_ant;                      %   Memory allocation

%   Calculates the likelihood vector
%likelihood_vect = measurement_likelihood(yk,xk_pred,Pk_pred,C_matrices,R_matrices,dimProbVect(1));

%   Calculates the Eta vector defined in "Online Bayesian Estimation of 
%   Transition Probabilities for Markovian Jump Systems"
Eta_vector = prob_vector_ant*( 1/( (prob_vector_ant')*TPM_ant*likelihood_vect ) );

%   Calculates the gij values and stores them into a matrix

for i=1:dimProbVect(1)
    for j=1:dimProbVect(1)      
        g_matrix(i,j) = 1 + Eta_vector(i,1)*(likelihood_vect(j,1)-TPM_ant(i,:)*likelihood_vect);        
    end
end

%Updates all the alpha parameters from the Dirichlet Distribution
for i=1:dimProbVect(1)    
    for j=1:dimProbVect(1)            
        %Denominator for the next operation
        sum_il = 0;
        for l=1:dimProbVect(1) 
            sum_il = sum_il+alpha_matrix_ant(i,l)*g_matrix(i,l);
        end        
        %Updates the alpha parameters from the Dirichlet Distribution
        alpha_matrix_pos(i,j) = alpha_matrix_ant(i,j)+(alpha_matrix_ant(i,j)*g_matrix(i,j))/sum_il;        
    end    
end

%Updates the TPM according to the measurement's likelihoods
for i=1:dimProbVect(1)    
    for j=1:dimProbVect(1)         
        TPM_pos(i,j) = (1/(iteration_number+inital_alpha_sum(i,1)))*alpha_matrix_pos(i,j);    
    end
end

%Takes the matrix transpose to make it organized in columns.
TPM_pos = TPM_pos';
alpha_matrix_pos = alpha_matrix_pos';




% %Gives the measurement's likelihood with respect to each possible mode. It
% %is assumed that the measurement equation is linear, yk = C*xk+w, and that
% %the measurement noise is Gaussian with zero mean.
% function [likelihood] = measurement_likelihood(yk,xk_pred,Pk_pred,C_matrices,R_matrices,dimProbVector)
%           
%     %The likelihood vector should have the same dimension as the mode
%     %probability vector
%     likelihood = zeros(dimProbVector,1);
% 
%     %dimProbVector = vector's dimension
%     for i=1:dimProbVector
% 
%         inov = yk - C_matrices(:,:,i)*xk_pred(:,:,i);
%         
%         V = C_matrices(:,:,i)*Pk_pred(:,:,i)*(C_matrices(:,:,i)')+R_matrices(:,:,i);
%         
%         %Auxiliary variable
%         C = 1/( ((2*pi)^(dimProbVector/2))*sqrt(det(V)) );
%         %Exponential's power
%         %power = (-0.5*(inov')*(V\inov));
%         power = 0.1;
%         
%         logLikelihood = power+log(C);
% 
%         % Returns the likelihood according to a multivariate Gaussian
%         % distribution
%         likelihood(i,1) = exp(logLikelihood);
% 
%     end





