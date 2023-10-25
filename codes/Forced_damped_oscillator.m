%% Haar Wavelet method to solve system of linear differential equations

                % KAUSHIK IYER
                
%           y(t) = x(t) + 0.3*x'(t), x(0) = 0
%           2.*y'(t) - 2 + 0.3*y(t) = 0.6(1-sin(t)),  y(0) = 0
             
%           the exact solution is 
%           x(t) = exp(-0.07*t)*((-278/109)*cos(0.703118*t) + (-110435000/38319931)*sin(0.703118*t)) + 2 + (200/109)*sin(t) + (60/109)*cos(t)      
clear;
clc;
addpath '/Users/kaushikiyer/Desktop/haar wavelets/codes/WHM'
% step 1
% collocation points
J    = input('Enter the level of decomposition:');  % level of decomposition
N    = 2^(J); % N = 2M                              % number of basis functions
j    = 1:N;                                         % index of grid points
t    = (j-0.5) ./ N;                                % grid points
alpha = 1;
M = 2; %mass
S = 1; %hook force constant
r = 0.3; % damping constant

%generating the Haar matrix and integral of Haar matrix of size (2M x 2M)
H = zeros(N,N);                                     % initialising the matrix with zeros
P = zeros(N,N);
for i = 1:N
    H(i,:) = haar_matrix(t,i,N); 
    P(i,:) = integral_of_H(t,i,N,alpha);
end
%disp(P);
D_alpha = eye(N)/(P*H');
a2 = zeros(1,N);
a1 = zeros(1,N);

% finding the filter coefficients
E  = ones(1,length(t)); 

a1 = fsolve(@(a1) fun(t,E,H,P,D_alpha,a1), zeros(1,N)) ;
% a2 = (a1*H)/P; 
% constructing the approximated solution 

approx_sol_y1 = a1*H;
approx_sol_y2 = (approx_sol_y1/H)*(P*H')*H;

% true solution

exact_sol_y2 = exp(-0.07.*t).*((-278/109).*cos(0.703118.*t) + (-110435000/38319931).*sin(0.703118.*t)) + 2.*E + (200/109).*sin(t) + (60/109).*cos(t);


% error
error = zeros(length(t),1);
error = abs(exact_sol_y2 - approx_sol_y2);

% Euclidean norm of the error
 
error_norm = norm(error,2);



 %% Plot graphics 
% fig:01
figure('color','w')
plot(t,approx_sol_y2,'g',t,exact_sol_y2,'rs')
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
function F = fun(t,E,H,P,D_alpha,a1)
F = (2.*a1*D_alpha*H + a1*P + 0.3*a1*H - 2.*(E-sin(t)));
end


