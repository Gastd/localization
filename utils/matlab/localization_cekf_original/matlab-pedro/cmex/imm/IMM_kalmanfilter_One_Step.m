%   Interacting Multiple Model  -   IMM
%

%   Linear model for each mode:
%
%       xk+1 = Ak(Mk)*xk + Bk(Mk)*uk + Dk(Mk)*wk
%       yk+1 = Ck+1(Mk)*xk+1 + vk+1
%
%       wk~N(0,Qk(Mk))  vk~N(0,Rk(Mk))  x0~N(x0_est,P0)
%

%   For this estimator, the Markov Chain's TPM matrix is such that
%   TPM(i,j)=P{xk+1=i|xk=j} -> The columns add to one!


function [xk_est_IMM,Pk_est_IMM,p_est_IMM,state_vectors_corr,cov_matrices_corr,state_vectors_pred,cov_matrices_pred] = IMM_kalmanfilter_One_Step(TPM,A_matrices,B_matrices,C_matrices,D_matrices,Q_matrices,R_matrices,u,y,state_vectors_ant,cov_matrices_ant,mode_prob_vector_ant)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% IMM Estimation

%Dimension of the initial state vector
dimState = size(state_vectors_ant);
%Dimension of the mode probability vector
dimMode = size(mode_prob_vector_ant);

%Memory allocation
xk_predictions = zeros(dimState(1),1,dimMode(1));
xk_corrections = zeros(dimState(1),1,dimMode(1));
Pk_predictions = zeros(dimState(1),dimState(1),dimMode(1));
Pk_corrections = zeros(dimState(1),dimState(1),dimMode(1));

%Updates the mode probability using the Markov TPM
[mode_prob_vector_pred] = markov_prob_update(mode_prob_vector_ant,TPM);

%Mixes the estimates
[mixed_state_vectors,mixed_cov_matrices] = mix_estimates(state_vectors_ant,cov_matrices_ant,mode_prob_vector_ant,mode_prob_vector_pred,TPM);

%Runs a Kalman filter matched to each mode

%dimMatricesA(3) = number of models
for i=1:dimMode(1)
    [xk_predictions(:,:,i), Pk_predictions(:,:,i), xk_corrections(:,:,i), Pk_corrections(:,:,i)] = kalmanFilter(A_matrices(:,:,i),B_matrices(:,:,i),C_matrices(:,:,i),D_matrices(:,:,i),Q_matrices(:,:,i),R_matrices(:,:,i),mixed_cov_matrices(:,:,i),mixed_state_vectors(:,:,i),u,y);
end

%Updates the mode probability
[updated_mode_prob_vector] = mode_prob_update(mode_prob_vector_pred,y,xk_predictions,Pk_predictions,C_matrices,R_matrices);

%Outputs the system's estimate
[xk_est_IMM,Pk_est_IMM] = imm_output(xk_corrections,Pk_corrections,updated_mode_prob_vector);
p_est_IMM = updated_mode_prob_vector;

%Outputs the state vectors for the different Kalman Filters
state_vectors_corr = xk_corrections;
cov_matrices_corr = Pk_corrections;

state_vectors_pred = xk_predictions;
cov_matrices_pred = Pk_predictions;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    
    


