%Updates all the mode probabilities for the Markov Chain. The Transition
%Probability Matrix (TPM) should be written column to row.
%
%   Mi,j = p(mk=i|mk-1=j)
%
function [pred_mode_prob_vector] = markov_prob_update(old_mode_prob_vector,Markov_matrix)
    
    pred_mode_prob_vector = Markov_matrix*old_mode_prob_vector;