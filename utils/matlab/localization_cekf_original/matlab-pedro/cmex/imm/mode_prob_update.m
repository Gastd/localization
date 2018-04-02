%Updates the mode probabilities according to the experimental measurement's 
%probability distribution, which is a Gaussian with parameters given by the
%prediction step from the multiple estimators.
function [updated_mode_prob_vector, likelihood_vector] = mode_prob_update(pred_mode_prob_vector,yk,yk_predictions,P_inov)
    
    %Retrieves the vector's dimension
    dimProbVector = size(pred_mode_prob_vector);      
    %Retrieves the output's dimension
    %dimOutput = size(yk);
    
    %Memory allocation
    updated_mode_prob_vector = zeros(dimProbVector);       
    likelihood_vector = zeros(dimProbVector);    
    %V = zeros(dimOutput(1),dimOutput(1),dimProbVector(1));
    
    debug2 = zeros(dimProbVector(1),1);
    
    power = zeros(dimProbVector(1),1);
    sqrtDetV = zeros(dimProbVector(1),1);
    
    %dimProbVector(1) = vector's dimension
    for i=1:dimProbVector(1)
        
        %Innovation
        inov = yk - yk_predictions(:,:,i);
        
        %Innovation's covariance matrix
        %V(:,:,i) = P_inov(:,:,i);                    
        V = P_inov(:,:,i);    
        
        %Power
        power(i,1) = -0.5*(inov')*(V\inov);       
                      
        debug2(i,1) = det(V);
                
        %Matrix's determinant square root
        sqrtDetV(i,1) = sqrt(det(V));
        
        %Measurement's likelihood
        C = 1/( ((2*pi)^(dimProbVector(1)/2))*sqrtDetV(i,1) );
        likelihood_vector(i,1) = C*exp(power(i,1));       
    end
    
    %likelihood_vector = [0.1;0.1];
    
    %Retrieves the minimium exponent
    %maxPower = max(power);    
    %maxSqrt = max(sqrtDetV);
    
    %These operations do not change the probabilities and avoid some
    %numerical issues related to very high (or very low) powers for the
    %exponential function and for the matrices' determinants.
    %power = power - ones(dimProbVector(1),1)*maxPower;
    %sqrtDetV = sqrtDetV*(1/maxSqrt);
    
    %Calculates the likelihood for all modes
    for i=1:dimProbVector(1)           
        updated_mode_prob_vector(i,1) = (pred_mode_prob_vector(i,1)/sqrtDetV(i,1))*exp(power(i,1));
    end
    
    %Calculates the probability sum
    prob_sum = sum(updated_mode_prob_vector);

    %Normalizes the vector
    updated_mode_prob_vector = updated_mode_prob_vector/prob_sum;
    
%     %This workaround solves the problem of very low likelihoods, where the
%     %exponential term can be 0. In this cases, the probability vector is
%     %reset.
%     if(prob_sum==0)
%         updated_mode_prob_vector = ones(dimProbVector)*1/(dimProbVector(1));
%     else
%         %Normalizes the vector
%         updated_mode_prob_vector = updated_mode_prob_vector/prob_sum;
%     end
    
    