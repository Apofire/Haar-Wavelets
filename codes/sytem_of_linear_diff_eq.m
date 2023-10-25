%% Haar Wavelet method to solve system of linear differential equations

                % KAUSHIK IYER
                
%           y1'(t) = y3(t) - cos(t), y1(0) = 1
%           y2'(t) = y3(t) - exp(t), y2(0) = 0
%           y3'(t) = y1(t) - y2(t),  y3(0) = 2
             
%           the exact solution is 
%        y1(t) = exp(t)     y2(t) = sin(t)      y3(t) = exp(t) + cos(t)
clear all;
addpath '/Users/kaushikiyer/Desktop/haar wavelets/codes/WHM'
% step 1
% collocation points
J    = input('Enter the level of decomposition:');  % level of decomposition
N    = 2^(J+1); % N = 2M                            % number of basis functions
j    = 1:N;                                         % index of grid points
t    = (2.*j - 1) ./ (2.*N);                                % grid points
alpha = 1;


%generating the Haar matrix and integral of Haar matrix of size (2M x 2M)
H = zeros(N,N);                                     % initialising the matrix with zeros
P = zeros(N,N);
for i = 1:N
    H(i,:) = haar_matrix(t,i,N); 
    P(i,:) = integral_of_H(t,i,N,alpha);
end
%disp(H);

% finding the filter coefficients
a1 = zeros(1,length(t));                            % coefficient function matrix at the collocation points
a2 = zeros(1,length(t));
E  = ones(1,length(t)); 

%a3 = fsolve(@(a3) fun(F_t,G_t,E,H,P,a3),zeros(1,length(t)));
% using the initial conditions to find the first estimate of the function
% and using it to find the filter coefficients.

a1 = (exp(t) - cos(t))/H;

% constructing the approximated solution 
approx_sol_y1_y2 = a1*P + E;

approx_sol_y3 = (approx_sol_y1_y2/H)*(P*H')*H + 2.*E;

a2 = ((approx_sol_y3) - exp(t))/H;
a3 = (approx_sol_y3 - cos(t))/H;

approx_sol_y2 = a2*P;
approx_sol_y1 = a3*P + E;

% true solution
exact_sol_y1 = exp(t);
exact_sol_y2 = sin(t);
exact_sol_y3 = exp(t) + cos(t);

% error
error = zeros(length(t),3);
error(:,1) = abs(exact_sol_y1 - approx_sol_y1);
error(:,2) = abs(exact_sol_y2 - approx_sol_y2);
error(:,3) = abs(exact_sol_y3 - approx_sol_y3);

% Euclidean norm of the error
error_norm = zeros(1,3);
error_norm(1,1) = norm(error(:,1),2); 
error_norm(1,2) = norm(error(:,2),2);
error_norm(1,3) = norm(error(:,3),2);


 %% Plot graphics 
% fig:01
figure('color','w')
plot(t,approx_sol_y1,'g',t,exact_sol_y1,'rs')
hold on
plot(t,approx_sol_y2,'g',t,exact_sol_y2,'rs')
hold on
plot(t,approx_sol_y3,'g',t,exact_sol_y3,'rs')
xlabel('$t$'); ylabel('$u(t)$');
title(['J = ' num2str(J) ', ' '2M = ' num2str(N)])
legend('Wavelet','Exact')


% fig:02
figure('color','w')
plot(t,error,'r.-')
xlabel('$t$'); ylabel('Absolute Error');
title('Absolute Error: $\max|y_{numeric} - y_{analytic}|$')





%% Functions

function y = haar_matrix(t,i,N)
% Function to generate the haar function for the given interval for a
% given index 'i'

% inputs
% x - the length of the vector containing the collocation points
% i - index of the haar function
% J - Maximum resolution
% outputs
% y = column vector containing the haar function
i = i-1;
if i ~= 0
   j = floor(log2(i));
   k = i - 2^j + 1; 
end

y = zeros([length(t) 1]);

if i == 0    
    for i = 1:length(t)
        if (0 <= t(i) && (t(i) < 1))
            y(i) = (1/sqrt(N))*1;
        else
            y(i) = 0;
        end
    end
else
    psi1 = (k -1) / (2^j);
    psi2 = (k - 0.5) / (2^j);
    psi3 = (k) / (2^j);
    for i = 1:length(t)
        if (psi1 <= t(i) && (t(i) < psi2))
            y(i) = (1/sqrt(N))*round(2^(j/2)*1,4);
        elseif (psi2 <= t(i) && (t(i) < psi3))
            y(i) = (1/sqrt(N))*round(2^(j/2)*(-1),4);
        else
            y(i) = 0;
        end
    end
end


end


function y = integral_of_H(t,i,N,alpha)
% Function to generate the integral of the haar function for the given interval for a
% given index 'i'

% inputs
% x - the length of the vector containing the collocation points
% i - index of the haar function
% J - Maximum resolution
% outputs
% y = column vector containing the integral of the haar function
i = i - 1; 
if i ~= 0
   j = floor(log2(i));
   k = i - 2^j + 1; 
end
    
y = zeros([1 length(t)]);

if i == 0    
    for i = 1:length(t)
        if (0 <= t(i) && (t(i) < 1))
            y(i) = (1/sqrt(N))*(1/(gamma(alpha+1)))*(t(i)^alpha);
        else
            y(i) = 0;
        end
    end
else
    psi1 = (k-1) / (2^j);
    psi2 = (k - 0.5) / (2^j);
    psi3 = (k) / (2^j);
    for i = 1:length(t)
        if (0 <= t(i) && (t(i) < psi1))
            y(i) = 0;
        elseif (psi1 <= t(i) && (t(i) < psi2))
            y(i) = (1/sqrt(N))*((2^(j/2))*((1/gamma(alpha+1))*(t(i) - psi1)^alpha));
        elseif (psi2 <= t(i) && (t(i) < psi3))
            y(i) = (1/sqrt(N))*((2^(j/2))*((1/gamma(alpha+1))*(t(i) - psi1)^alpha - (2/gamma(alpha+1))*(t(i) - psi2)^alpha)) ;
        elseif (psi3 <= t(i) && (t(i) < 1))
            y(i) = (1/sqrt(N))*((2^(j/2))*((1/gamma(alpha+1))*(t(i) - psi1)^alpha - (2/gamma(alpha+1))*(t(i) - psi2)^alpha + (1/gamma(alpha+1))*(t(i) - psi3)^alpha));   
        end
    end
end
end

function F = fun(F_t,G_t,E,H,P,a3)
F = (a3*(H) - (F_t - G_t)*(H\P) - E);
end


