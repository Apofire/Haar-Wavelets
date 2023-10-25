%% Haar Wavelet method to solve Burger-Huxley equation

                % KAUSHIK IYER
                

addpath '/Users/kaushikiyer/Desktop/haar wavelets/codes/WHM'

% collocation points
J    = input('Enter the level of decomposition:');  % level of decomposition
N    = 2^(J+1); % N = 2M                            % number of basis functions
j    = 1:N;                                         % index of grid points
x    = (2*j-1) ./ (2*N);                            % grid points
t   = 0;                                            % initial value
                                                    % change of variable to [0,0.5]
delta_t = 0.1;                                  % time increments 
iter = input('Enter number of iterations of time:');
alpha = input('Enter value of the constant alpha:');
beta = input('Enter value of the constant beta:');
gamma = input('Enter value of the constant gamma:');
delta = input('Enter value of the constant delta:');

a1 = ((-alpha*delta + delta*sqrt(alpha^2 + 4*beta*(1 + delta)))/(4*(1 + delta)))*gamma;
a2 = alpha*gamma/(1 + delta) - ((1 + delta - gamma)*(-alpha + sqrt(alpha^2 + 4*beta*(1 + delta))))/(2*(1 + delta));

% generating the Haar matrix and its integral matrices of size (2M x 2M)
H = zeros(N,N);   % initialising the matrix with zeros
Q = zeros(N,N);   % initialising the matrix with zeros
P = zeros(N,N);
Q_1 = zeros(N,N); % value of Q(1)

for i = 1:N
    H(i,:) = haar_matrix(x,i,J); 
    Q(i,:) = integral_of_P(x,i,J);
    P(i,:) = integral_of_H(x,i,J);
    Q_1(i,:) = q_1(x,i,J,1,0);
end
%disp(H);

E = ones(1,N);

% defining initial and boundary conditions
f_x = (gamma/2.*E + gamma/2*tanh(a1.*x)).^(1/delta);    
f_x_d1 = (a1/delta)*(gamma/2)^(1/delta).*(sech(a1.*x).^2).*(E + tanh(a1.*x).^(1/delta) -1);                       % first derivative of f_x
f_x_d2 = -(a1^2/delta^2)*((gamma/2)^(1/delta)).*(sech(a1.*x).^2).*((E + tanh(a1.*x)).^(1/delta -1)).*(2*delta.*(E + tanh(a1.*x).*tanh(a1.*x) + (delta-1).*(sech(a1.*x).^2)));
% g0_t = (gamma/2 + gamma/2*tanh(-(a1*a2)*t))^(1/delta); 
% g0_d1_t  = (a1*a2*gamma*(tanh(a1*a2*t)^2 - 1)*(gamma/2 - (gamma*tanh(a1*a2*t))/2)^(1/delta - 1))/(2*delta);   
% 
% g1_t = (gamma/2 + gamma/2*tanh(a1(1 - a2*t)))^(1/delta);                 
% g1_d1_t = (a1*a2*gamma*(tanh(a1*(a2*t - 1))^2 - 1)*(gamma/2 - (gamma*tanh(a1*(a2*t - 1)))/2)^(1/delta - 1))/(2*delta);

a = zeros(iter,N);

u_ts = f_x;
u_d1_ts = f_x_d1;
u_d2_ts = f_x_d2;




approx_sol = zeros(iter,N);

% finding the filter coefficients
for ts = 1:iter
    
    % updating the equations
    if ts ~=1
        u_d2_ts = delta_t.*(a(ts,:)*H) + u_d2_ts;
        u_d1_ts = delta_t.*a(ts,:)*(P - E.*Q_1) + u_d1_ts + (g1_ts1 - g1_ts + g0_ts - g0_ts1).*E;
    end
    t = t + delta_t;            % time increment 
    
    g0_ts1     = (gamma/2 + gamma/2*tanh(-(a1*a2)*t))^(1/delta); 
    g0_ts      = (gamma/2 + gamma/2*tanh(-(a1*a2)*(t - delta_t)))^(1/delta);
    g0_d1_ts1  = (a1*a2*gamma*(sech(a1*a2*t)^2)*(gamma/2 - (gamma*tanh(a1*a2*t))/2)^(1/delta - 1))/(2*delta);   
    
    g1_ts1     = (gamma/2 + gamma/2*tanh(a1*(1 - a2*t)))^(1/delta); 
    g1_ts      = (gamma/2 + gamma/2*tanh(a1*(1 - a2*(t - delta_t))))^(1/delta);            
    g1_d1_ts1  = (a1*a2*((gamma/2)^(1/delta))*((sech(a1 - a1*a2*t))^2)*(tanh(a1 - a1*a2*t)^(1/delta - 1)))/delta;

    % finding filter coefficients
    %a = fsolve(@(a) myfun(a,x,alpha,beta,gamma,delta,u_ts,g0_d1_ts1,g1_d1_ts1,E,Q,Q_1,u_d2_ts,u_d1_ts), 0.1.*ones(1,N));
    a(ts,:) = (u_d2_ts - alpha.*((u_ts.^(delta)).*u_d1_ts) + beta.*u_ts.*(E - u_ts.^delta).*(u_ts.^delta - gamma.*E) - x.*(g0_d1_ts1 - g1_d1_ts1) - g0_d1_ts1.*E)/(Q - x.*Q_1);
   
    % constructing the approximated solution 
    approx_sol(ts,:) = (delta_t).*(a(ts,:)*(Q - x.*Q_1)) + u_ts + (g0_ts1 - g0_ts).*E - x.*(-g1_ts1 + g1_ts - g0_ts + g0_ts1); 

    u_ts = approx_sol(ts,:); % saving the solution of the previous time index 

    
   
    
 
end
t = 0;
% true solution
exact_sol = zeros(iter,N);
for i = 1:iter
    t = t + delta_t;    
    exact_sol(i,:) = ((gamma/2).*E + (gamma/2).*(tanh(a1.*x - (a1*a2*t).*E))).^(1/delta);
     
end

% error
error = zeros(iter,N);

for i = 1:iter
    error(i,:) = abs(exact_sol(i,:) - approx_sol(i,:));
    %error_norm = norm(error,2); % Euclidean norm of the error
end


%% Plot graphics
set(0,'defaulttextinterpreter','latex')
set(0,'defaultaxesfontname','times')
set(0,'defaultaxesfontsize',12)

figure(1)
surf(approx_sol)
hold on
surf(exact_sol)
xlabel('$x$'); ylabel('$t$'); zlabel('$y(x,t)$');
title(['J = ' num2str(J) ', ' '2M = ' num2str(N)])
legend('Wavelet','Exact')


figure(2)
surf(error)
xlabel('$x/32$'); ylabel('$t/10$'); zlabel('Absolute Error');
title('Absolute Error: $\max|y_{numeric} - y_{analytic}|$')




%% Functions

function y = haar_matrix(t,i,J)
% Function to generate the haar function for the given interval for a
% given index 'i'

% inputs
% x - the length of the vector containing the collocation points
% i - index of the haar function
% J - Maximum resolution
% outputs
% y = column vector containing the haar function

if i == 1
    m = 0;
    k = 0;
    j = 0;
else
    IMat = zeros([J+1 2^J]);
    IMask = IMat;

    ind_s = 1;
    for ind_j = 0:J    
        for ind_i = 0:2^ind_j-1
            ind_s = ind_s + 1;
            IMask(ind_j+1,ind_i+1) = ind_s;     
            IMat(ind_j+1,ind_i+1) = ind_i+ind_j+1;
        end % for i
    end % for j

    [ind_j, ind_i] = find(IMask == i);
    m = 2^(ind_j - 1);
    k = ind_i - 1;
    j = ind_j -1;
end

y = zeros([length(t) 1]);

if i == 1    
    for i = 1:length(t)
        if (0 <= t(i) && (t(i) < 1))
            y(i) = 1;
        else
            y(i) = 0;
        end
    end
else
    alpha = k / m;
    beta = (k + 0.5) / m;
    gamma = (k + 1) / m;
    for i = 1:length(t)
        if (alpha <= t(i) && (t(i) < beta))
            y(i) = 1;
        elseif (beta <= t(i) && (t(i) < gamma))
            y(i) = -1;
        else
            y(i) = 0;
        end
    end
end


end


function y = integral_of_H(t,i,J)
% Function to generate the integral of the haar function for the given interval for a
% given index 'i'

% inputs
% x - the length of the vector containing the collocation points
% i - index of the haar function
% J - Maximum resolution
% outputs
% y = column vector containing the integral of the haar function

if i == 1
    m = 0;
    k = 0;
else
    IMat = zeros([J+1 2^J]);
    IMask = IMat;

    ind_s = 1;
    for ind_j = 0:J    
        for ind_i = 0:2^ind_j-1
            ind_s = ind_s + 1;
            IMask(ind_j+1,ind_i+1) = ind_s;     
            IMat(ind_j+1,ind_i+1) = ind_i+ind_j+1;
        end % for i
    end % for j

    [m, k] = find(IMask == i);
    m = 2^(m - 1);
    k = k - 1;
end


y = zeros([1 length(t)]);

if i == 1    
    for i = 1:length(t)
        if (0 <= t(i) && (t(i) < 1))
            y(i) = t(i);
        else
            y(i) = 0;
        end
    end
else
    alpha = k / m;
    beta = (k + 0.5) / m;
    gamma = (k + 1) / m;
    for i = 1:length(t)
        if (alpha <= t(i) && (t(i) < beta))
            y(i) = t(i) - alpha;
        elseif (beta <= t(i) && (t(i) < gamma))
            y(i) = gamma - t(i);
        else
            y(i) = 0;
        end
    end
end
end


function y = integral_of_P(t,i,J)
% Function to generate the integral of the P function for the given interval for a
% given index 'i'

% inputs
% x - the length of the vector containing the collocation points
% i - index of the haar function
% J - Maximum resolution
% outputs
% y = column vector containing the integral of the P function

if i == 1
    m = 0;
    k = 0;
else
    IMat = zeros([J+1 2^J]);
    IMask = IMat;

    ind_s = 1;
    for ind_j = 0:J    
        for ind_i = 0:2^ind_j-1
            ind_s = ind_s + 1;
            IMask(ind_j+1,ind_i+1) = ind_s;     
            IMat(ind_j+1,ind_i+1) = ind_i+ind_j+1;
        end % for i
    end % for j

    [m, k] = find(IMask == i);
    m = 2^(m - 1);
    k = k - 1;
end


y = zeros([1 length(t)]);

if i == 1    
    for i = 1:length(t)
        if (0 <= t(i) && (t(i) < 1))
            y(i) = 0.5 * t(i) * t(i);
        else
            y(i) = 0;
        end
    end
else
    alpha = k / m;
    beta = (k + 0.5) / m;
    gamma = (k + 1) / m;
    for i = 1:length(t)
        if (alpha <= t(i) && (t(i) < beta))
            y(i) = vpa(0.5 * (t(i) - alpha).^2);
        elseif (beta <= t(i) && (t(i) < gamma))
            y(i) = vpa(1 / (4*m^2) - 0.5 * (gamma - t(i)).^2);
        elseif (gamma <= t(i) && (t(i) < 1))
            y(i) = vpa(1 / (4*m^2));
        else
            y(i) = 0;
        end
    end
end
end


function y = q_1(t,i,J,A,B)
% Function to generate the Q_1 function for the given interval for a
% given index 'i'
% It is defined as            [0.5(B-A)^2 , i=1  
%                     Q_i(1) =[
%                             [ 0.25*(B-A)^2/(m^2) , i>1

% inputs
% x - the length of the vector containing the collocation points
% i - index of the haar function
% J - Maximum resolution
% outputs
% y = column vector containing the Q_1 function

if i == 1
    m = 0;
    k = 0;
else
    IMat = zeros([J+1 2^J]);
    IMask = IMat;

    ind_s = 1;
    for ind_j = 0:J    
        for ind_i = 0:2^ind_j-1
            ind_s = ind_s + 1;
            IMask(ind_j+1,ind_i+1) = ind_s;     
            IMat(ind_j+1,ind_i+1) = ind_i+ind_j+1;
        end % for i
    end % for j

    [m, k] = find(IMask == i);
    m = 2^(m - 1);
    k = k - 1;
end

y = zeros([1 length(t)]);

    if i == 1       
      for i = 1:length(t)
       y(i) = 0.5*(B-A)^2;     
    end
else
    for i = 1:length(t)
         y(i) = (B-A)^2 / (4*m^2);
    end
    end
end


function F = myfun(a,x,alpha,beta,gamma,delta,u_ts,g0_d1_ts1,g1_d1_ts1,E,Q,Q_1,u_d2_ts,u_d1_ts)
F = a*(Q - x.*Q_1) - u_d2_ts + alpha.*((u_ts.^(delta)).*u_d1_ts) - beta.*u_ts.*(E - u_ts.^delta).*(u_ts.^delta - gamma.*E) + x.*(g0_d1_ts1 - g1_d1_ts1) + g0_d1_ts1.*E;

end





