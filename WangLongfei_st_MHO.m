%% 
%MILSTRONG Test strong convergence of Milstein: vectorized
%
% Solves dX = r*X*(K-X) dt + beta*X dW, X(0) = Xzero,
% where r = 2, K= 1, beta = 1 and Xzero = 0.5.
%
% Discretized Brownian path over [0,1] has dt = 2^(-11).
% Milstein uses timesteps 128*dt, 64*dt, 32*dt, 16*dt (also dt for reference).
%
% Examines strong convergence at T=1: E | X_L - X(T) |.
% Code is vectorized: all paths computed simultaneously.
rand('state',100)
r = 2; K = 1; beta = 0.25; Xzero = 0.5; % problem parameters
T = 1; N = 2^(11); dt = T/N; %
M = 500; % number of paths sampled
R = [1; 16; 32; 64; 128]; % Milstein stepsizes are R*dt
dW = sqrt(dt)*randn(M,N); % Brownian increments
Xmil = zeros(M,5); % preallocate array
for p = 1:5
Dt = R(p)*dt; L = N/R(p); % L timesteps of size Dt = R dt
Xtemp = Xzero*ones(M,1);
for j = 1:L
Winc = sum(dW(:,R(p)*(j-1)+1:R(p)*j),2);
Xtemp = Xtemp + Dt*r*Xtemp.*(K-Xtemp) + beta*Xtemp.*Winc ...
+ 0.5*beta^2*Xtemp.*(Winc.^2 - Dt);
end
Xmil(:,p) = Xtemp; % store Milstein solution at t =1
end
Xref = Xmil(:,1); % Reference solution
Xerr = abs(Xmil(:,2:5) - repmat(Xref,1,4)); % Error in each path
mean(Xerr); % Mean pathwise erorrs
Dtvals = dt*R(2:5); % Milstein timesteps used
subplot(224) % lower RH picture
loglog(Dtvals,mean(Xerr),'b*-'), hold on
loglog(Dtvals,Dtvals,'r--'), hold off % reference slope of 1
axis([1e-3 1e-1 1e-4 1])
xlabel('\Delta t')
ylabel('Sample average of | X(T) - X_L |')
title('milstrong.m','FontSize',10)
%%%% Least squares fit of error=C* Dt^q %%%%
A = [ones(4,1), log(Dtvals)]; rhs = log(mean(Xerr)');
sol = A\rhs; q = sol(2)
resid = norm(A*sol - rhs)
%% 


clear
clc
% initial parameters
l=0.2; r=0.033; v=1; %set the velocity
w1=v/r; w2=v/r;%two omega are same
ti=2^11;% times(how many pathes)
T=1; %total time
dt=T/ti; %time step
D=4; %noise coeffcient
P=T/dt; %how many points in one path
R = [1; 16; 32; 64; 128];
for p = 1:5
Dt = R(p)*dt; L = P/R(p); % L timesteps of size Dt = R dt

%Brownian increments
randn('state',400)
dw1=sqrt(Dt)*randn(ti,P);% two omega increase randomly and differently
dw2=sqrt(Dt)*randn(ti,P);
x = zeros(ti,L);      
y = zeros(ti,L);   
theta=zeros(ti,L);
%plot the distribution in cartesian
for i=1:ti %for the first path
    for j=2:L % from start time to end
        x(i,j) = x(i,j-1) + 0.5*r*(w1+w2)*cos(theta(i,j-1))*Dt + sqrt(D)*0.5*r*cos(theta(i,j-1))*(dw1(i,j-1)+dw2(i,j-1))-0.5*sqrt(D)*0.5*r*cos(theta(i,j-1))*sqrt(D)*0.5*r*sin(theta(i,j-1))*(dw1(i,j-1)^2+dw2(i,j-1)^2-Dt);
        y(i,j) = y(i,j-1) + 0.5*r*(w1+w2)*sin(theta(i,j-1))*Dt + sqrt(D)*0.5*r*sin(theta(i,j-1))*(dw1(i,j-1)+dw2(i,j-1))+0.5*sqrt(D)*0.5*r*sin(theta(i,j-1))*sqrt(D)*0.5*r*cos(theta(i,j-1))*(dw1(i,j-1)^2+dw2(i,j-1)^2-Dt);
        theta(i,j) = theta(i,j-1) + Dt*r*(w1-w2)/l + sqrt(D)*r*(dw1(i,j-1)-dw2(i,j-1))/l;
    end
%     plot(x(i,P),y(i,P),'b.')%plot every end point of evry path
%     hold on
end
Xmil(:,p) = x(:,end); % store Milstein solution at t =1
Ymil(:,p) = y(:,end);
end
Xref = Xmil(:,1); % Reference solution
Yref = Ymil(:,1);
Xerr=abs(Xmil(:,2:5) - repmat(Xref,1,4));
% err = sqrt((Xmil(:,2:5) - repmat(Xref,1,4)).^2+(Ymil(:,2:5)-repmat(Yref,1,4)).^2); % Error in each path
err = abs(Xmil(:,2:5) - repmat(Xref,1,4))+abs(Ymil(:,2:5)-repmat(Yref,1,4)); % Error in each path
Dtvals = dt*R(2:5); % Milstein timesteps used
figure;
loglog(Dtvals,mean(Xerr),'b*-'), hold on
% loglog(Dtvals,Dtvals,'r--'), hold off % reference slope of 1
% axis([1e-2 1e-1 1e-4 1])
xlabel('\Delta t')
ylabel('Sample average of | X(T) - X_L |')
title('Milstrong','FontSize',10)


