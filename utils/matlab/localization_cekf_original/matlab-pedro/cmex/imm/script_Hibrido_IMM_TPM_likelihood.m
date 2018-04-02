%   Simulation script for Hibrido_IMM_TPM.m
%
%   
%   This script simulates the tracking of an unidimensional state vector
%   (a water column's height, for example) using two sensors subject to 
%   heavy disturbances.
%
%
%
%   Author: Pedro Henrique R.Q.A. Santana


%% System model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%   Linear model for each mode:
%
%       xk+1 = Ak(Mk)*xk + Bk(Mk)*uk + Dk(Mk)*wk
%       yk+1 = Ck+1(Mk)*xk+1 + vk+1
%
%       wk~N(0,Qk(Mk))  vk~N(0,Rk(Mk))  x0~N(x0_est,P0)
%

clear all;

%Number of iterations
num_iter = 1000;

%Markov chain matrix
%
%   m1 = sensor s1 is operating normally
%   m2 = sensor s2 is operating normally
%   m3 = both sensors are operating normally
%   m3 = both sensors are faulty
%
%   Column sum = 1

TPM = [0.6 0.1 0.15 0.15; 0.1 0.6 0.15 0.15; 0.15 0.15 0.6 0.1; 0.15 0.15 0.1 0.6];

%System matrix
A = 1;      % The prediction step is the same for all the modes    

%Input matrix
B = 0;      % There is no input for all the modes.

%Output matrices
C1 = [1;0];
C2 = [0;1];
C3 = [1;1];
C4 = [0;0];

%Model noise matrix
D = 1;

%Model covariance matrix
Q = 10;

%Sensor's output variance (working)
r1 = 0.2;   %   Sensor 1
r2 = 0.5;     %   Sensor 2

%Sensor's output variance (fault)
f1 = 278;   %   Fault in Sensor 1
f2 = 278;   %   Fault in Sensor 2

%Output covariance matrices for each state
R1 = [r1 0;0 f2];
R2 = [f1 0;0 r2];
R3 = [r1 0;0 r2];
R4 = [f1 0;0 f2];

%Organization of all system's models into a data structure

A_matrices = zeros(1,1,4);
for i=1:4
    A_matrices(:,:,i) = A;
end

B_matrices = zeros(1,1,4);
for i=1:4
    B_matrices(:,:,i) = B;
end

C_matrices = zeros(2,1,4);
C_matrices(:,:,1) = C1;
C_matrices(:,:,2) = C2;
C_matrices(:,:,3) = C3;
C_matrices(:,:,4) = C4;

D_matrices = zeros(1,1,4);
for i=1:4
    D_matrices(:,:,i) = D;
end

Q_matrices = zeros(1,1,4);
for i=1:4
    Q_matrices(:,:,i) = Q;
end

R_matrices = zeros(2,2,4);
R_matrices(:,:,1) = R1;
R_matrices(:,:,2) = R2;
R_matrices(:,:,3) = R3;
R_matrices(:,:,4) = R4;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Simulation

%   Lets simulate an object that follows a random walk in a line from
%   [-50,50]. In order to keep track of this particle, two radars subject
%   to faulty measurements are used. The faults occur according to a Markov
%   Chain, whose real TPM is unknown.

%Initial conditions
x0 = 0;             %State%
P0 = 11;            %State variance%
p0 = [0;0;1;0];     %Mode%

%Random Walk parameters
rw_var = 3;     %   Random Walk's variance. This should be at least 3 times
                %   smaller than Q.

%Markov Chain Parameters
mode = 3;                       %   Initial mode (Both sensors working)
TPM_ant = ones(4,4)*0.25;       %   Initial TPM (Complete ignorance)
alpha_matrix_ant = ones(4,4);   %   Initial values for the Dirichlet Distributions' parameters


%Initial state
xk_ant_sim = x0;

%Simulated measurements, states and modes
y_sim = zeros(2,1,num_iter);
x_sim = zeros(num_iter,1);
m_sim = zeros(num_iter,1);


%Initial values
x_sim(1,1) = x0;
y_sim(:,:,1) = x0;
m_sim(1,1) = mode;

xk_est_IMM_ant = zeros(1,1,4);
xk_est_IMM_ant(:,:,1) = x0;
xk_est_IMM_ant(:,:,2) = x0;
xk_est_IMM_ant(:,:,3) = x0;
xk_est_IMM_ant(:,:,4) = x0;

Pk_est_IMM_ant = zeros(1,1,4);
Pk_est_IMM_ant(:,:,1) = P0;
Pk_est_IMM_ant(:,:,2) = P0;
Pk_est_IMM_ant(:,:,3) = P0;
Pk_est_IMM_ant(:,:,4) = P0;

p_est_IMM_ant = p0;

xk_kalman_corr = zeros(num_iter,1);
xk_kalman_corr(1,1) = x0;

Pk_kalman_corr = zeros(1,1,num_iter);
Pk_kalman_corr(:,:,1) = P0;

xk_kalman_ant = x0;
Pk_kalman_ant = P0;


%Data storage variables
xk_est_IMM = zeros(num_iter,1);
xk_est_IMM(1,1) = x0;

Pk_est_IMM = zeros(1,1,num_iter);
Pk_est_IMM(:,:,1) = P0;

p_est_IMM = zeros(4,1,num_iter);
p_est_IMM(:,:,1) = p0;

%   In this simulation, perfect mode measurement is supposed an in
%   "Estimation of Non-sationary Markov Chain Transition Models"
for k=2:num_iter
          
    %sprintf('Iteracao = %d\n',k)
    
    %Simulates the markov chain's transitions    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    transition_prob = rand(1);    
    for i=1:4        
        if transition_prob <= TPM(i,mode)
            mode = i;
            break;
        else
            transition_prob = transition_prob-TPM(i,mode);  
        end
    end    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Simulates a state evolution and a new measurement
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %System's evolution according to the random walk
    x_sim(k,1) = xk_ant_sim+randn(1)*rw_var;    
    %Measurement according to the mode
    y_sim(:,:,k) = C_matrices(:,:,mode)*x_sim(k,1) + chol(R_matrices(:,:,mode))*randn(2,1);        
    %Stores the mode    
    m_sim(k,1) = mode;    
    %Updated the previous state
    xk_ant_sim = x_sim(k,1);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %State and mode estimation according to the IMMKF
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %   IMMKF
    [xk_est_IMM(k,1),Pk_est_IMM(:,:,k),p_est_IMM(:,:,k),xk_est_IMM_corr,Pk_est_IMM_corr,xk_est_IMM_pred,Pk_est_IMM_pred] = IMM_kalmanfilter_One_Step(TPM_ant,A_matrices,B_matrices,C_matrices,D_matrices,Q_matrices,R_matrices,0,y_sim(:,:,k),xk_est_IMM_ant,Pk_est_IMM_ant,p_est_IMM_ant);
                        
    %   Likelihood Dirichlet Estimator
    [TPM_pos,alpha_matrix_pos] = likelihoodDirichletEstimator(y_sim(:,:,k),xk_est_IMM_pred,Pk_est_IMM_pred,C_matrices,R_matrices,TPM_ant,p_est_IMM_ant,alpha_matrix_ant,k-1,[4;4;4;4]);
    
    %Updates the variables for a new iteration
    TPM_ant = TPM_pos;
    alpha_matrix_ant = alpha_matrix_pos;
    
    xk_est_IMM_ant = xk_est_IMM_corr;
    Pk_est_IMM_ant = Pk_est_IMM_corr;
    p_est_IMM_ant = p_est_IMM(:,:,k);
    
           
    %Standard Kalman Filter
    [xk_kalman_pred, Pk_kalman_pred, xk_kalman_corr(k,1), Pk_kalman_corr(:,:,k)] = kalmanFilter(A,B,C3,D,Q,R3,Pk_kalman_ant,xk_kalman_ant,0,y_sim(:,:,k));
    
    xk_kalman_ant = xk_kalman_corr(k,1);
    Pk_kalman_ant = Pk_kalman_corr(:,:,k);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figures
close all;

%   Removes the singleton dimension
y_sim = squeeze(y_sim);
p_est_IMM = squeeze(p_est_IMM);

%   Horizontal axis
hor_axis = 1:1:num_iter;

%   IMM state estimation
figure;
plot(hor_axis,x_sim,'b',hor_axis,xk_est_IMM,'r');
legend('Real','IMM');
ylabel('x(k)');
xlabel('k');
title('Hybrid IMM estimator');

%   KF state estimation
figure;
plot(hor_axis,x_sim,'b',hor_axis,xk_kalman_corr,'r');
legend('Real','KF');
ylabel('x(k)');
xlabel('k');
title('Kalman Filter');

%   Measurements from Sensor 1
figure;
plot(hor_axis,x_sim,'b',hor_axis,y_sim(1,:),'r');
legend('Real','Sensor 1');
ylabel('x(k)');
xlabel('k');
title('Sensor 1');

%   Measurements from Sensor 2
figure;
plot(hor_axis,x_sim,'b',hor_axis,y_sim(2,:),'r');
legend('Real','Sensor 2');
ylabel('x(k)');
xlabel('k');
title('Sensor 2');

%Prints the real and estimated TPMs
TPM=TPM
TPM_pos = TPM_pos


% figure;
% plot(hor_axis,x_sim(2,:),'r',hor_axis,xk_est_IMM(2,:),'b');
% legend('Real','Estimation');
% ylabel('$$x_{2}(k)$$','Interpreter','latex','FontSize',25);
% xlabel('k','Interpreter','latex','FontSize',25);
% title('State estimation results','Interpreter','latex','FontSize',25);
% 
% %Mode
% figure;
% plot(hor_axis,m_est_IMM(1,:),'r',hor_axis,m_est_IMM(2,:),'b');
% legend('Mode 1','Mode 2');
% ylabel('$$p(k)$$','Interpreter','latex','FontSize',25);
% xlabel('k','Interpreter','latex','FontSize',25);
% title('State estimation results','Interpreter','latex','FontSize',25);
% 
% %Determines the most likely mode
% mode_est = zeros(num_iter,1);
% 
% for i=1:num_iter
%     if m_est_IMM(1,i)>0.5
%         mode_est(i,1) = 1;
%     else
%         mode_est(i,1) = 2;
%     end
% end
% 
% figure;
% plot(hor_axis,m_sim,'rs',hor_axis,mode_est,'bx');
% legend('Real','Estimation');
% ylabel('$$m(k)$$','Interpreter','latex','FontSize',25);
% xlabel('k','Interpreter','latex','FontSize',25);
% title('State estimation results','Interpreter','latex','FontSize',25);















