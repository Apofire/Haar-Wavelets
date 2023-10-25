%% Haar Wavelet method to solve fractional partial differential equations

                % KAUSHIK IYER
                
%  D^(α)u(x,t) + D^(β)*u(x,t) = f(x,t)    where D^(α) is wrt x and D^(β) is wrt t       
% with the initial condition u(0) = 0    alpha = 1/3 , beta = 1/2      
% f(x,t) = 9*x^2*(t^5/3)/(5*gamma(2/3)) + 8*(x^3/2)*t^2/(3*gamma(1/2))
%           the exact solution is 
%        u(x,t) = x^2*t^2
clear;
clc;
addpath '/Users/kaushikiyer/Desktop/haar wavelets/codes/WHM'
% step 1
% collocation points
J    = input('Enter the level of decomposition:');  % level of decomposition
N    = 2^(J); % N = 2M                            % number of basis functions
j    = 1:N;                                         % index of grid points
x    = (j-0.5) ./ N;                                % grid points
alpha = input('Enter the value of alpha:');         % order of differential for space
beta = input('Enter the value of beta:');           % order of differential for time
t    = (j-0.5)./N;

%generating the Haar matrix of size (2M x 2M)
H = zeros(N,N);                                     % initialising the matrix with zeros
for i = 1:N
    H(i,:) = haar_matrix(x,i,N);
    P_alpha_initial(i,:) = integral_of_H(x,i,N,alpha);
    P_beta_initial(i,:) = integral_of_H(x,i,N,beta);
end
% disp(P);
P_alpha = P_alpha_initial*H';
P_beta  = P_beta_initial*H';
%disp(P_alpha);

% generating the operational matrix of fractioanl differentiation (NxN)
D_alpha = round(eye(N)/P_alpha,10);
D_beta  = round(eye(N)/P_beta,10);
%disp(D_alpha);

approx_sol = zeros(length(t),length(x));
exact_sol  = zeros(length(t),length(x)); 
error      = zeros(length(t),length(x));


f = zeros(length(t),length(x));

% finding the filter coefficients
for i = 1:length(x)
    for j = 1:length(t)
        f(i,j) = (9*(x(i)^2)*(t(j)^5/3))/(5*gamma(2/3)) + (8*(x(i)^3/2)*t(j)^2)/(3*gamma(1/2));    % function values at the collocation points
        
        % true solution
        exact_sol(i,j) = (x(i)^2)*(t(j)^2);
    end
end

c1 = H*f*H';
     c  = fsolve(@(c) myfun(c,H,c1,D_alpha,D_beta), zeros(N,N));
%     c1 = H*exact_sol*H';

    % constructing the approximated solution 
    approx_sol = H'*c*H;
    
   

    % error
    error = abs(exact_sol - approx_sol);
    
   


%% Plot graphics
set(0,'defaulttextinterpreter','latex')
set(0,'defaultaxesfontname','times')
set(0,'defaultaxesfontsize',12)

% fig:01
figure('color','w')
surf(approx_sol)
hold on
surf(exact_sol)
xlabel('$x$'); ylabel('$t$');zlabel('$u(x,t)$');
title(['J = ' num2str(J) ', ' '2M = ' num2str(N)])
legend('Wavelet','Exact')



% fig:02
figure('color','w')
surf(error)
xlabel('$x$');ylabel('$t$'); zlabel('Absolute Error');
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




function y = myfun(c,H,c1,D_alpha,D_beta)
y = D_alpha'*c + c*D_beta - c1;
end
