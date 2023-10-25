%% Haar Wavelet method to solve system of linear differential equations

                % KAUSHIK IYER
                
%           y1'(t) = y1(t) - 3y2(t), y1(0) = 0
%           y2'(t) = -3y2(t) + 5y3(t) , y2(0) = 0
%           y3'(t) = -5y3(t) -5*z0,  y3(0) = z0 = 1
             
%           the exact solution is 
%        y1(t) = (15*z0/8)*(exp(-5t) -2*exp(-3t) + exp(-t))  
%        y2(t) = (5*z0/2)*(-exp(-5t) + exp(-3t))
%        y3(t) = z0*(exp(-5t))

addpath '/Users/kaushikiyer/Desktop/haar wavelets/codes/WHM'
% step 1
% collocation points
J    = input('Enter the level of decomposition:');  % level of decomposition
N    = 2^(J); % N = 2M                            % number of basis functions
j    = 1:N;                                         % index of grid points
t    = (j-0.5) ./ N;                                % grid points
z0 = input('Enter initial unit of living tree biomass:');
alpha = 1;

%generating the Haar matrix and integral of Haar matrix of size (2M x 2M)
H = zeros(N,N);                                     % initialising the matrix with zeros
P = zeros(N,N);
for i = 1:N
    H(i,:) = haar_matrix(t,i,N); 
    P(i,:) = integral_of_H(t,i,N,alpha);
end
%disp(H);
P_alpha = P*H';
D_alpha = eye(N)/P_alpha;

% finding the filter coefficients
a1 = zeros(1,length(t));                            % coefficient function matrix at the collocation points
a2 = zeros(1,length(t));
a3 = zeros(1,length(t));
E  = ones(1,length(t)); 

a3 = (-5.*E)/(D_alpha*H + 5.*(H));
approx_sol_y3 = a3*H;

a2 = (5.*(a3*H+E))/(D_alpha*H + 3.*H);
approx_sol_y2 = a2*H;


a1 = (-3.*(a2*H))/(D_alpha*H - H);
approx_sol_y1 = a1*H;




% constructing the approximated solution 


% true solution
exact_sol_y1 = (15*z0/8).*(exp(-5.*t) -2*exp(-3.*t) + exp(-t));
exact_sol_y2 = (5*z0/2).*(-exp(-5.*t) + exp(-3.*t));
exact_sol_y3 = z0.*(exp(-5.*t));

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



