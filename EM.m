%EM Euler-Maruyama method on linear SDE
%
% SDE is dX = lambda*X dt + mu*X dW, X(0) = Xzero,
% where lambda = 2, mu = 1 and Xzero = 1.
%
% Discretized Brownian path over [0,1] has dt = 2^(-8).
% Euler-Maruyama uses timestep R*dt.
randn('state',100)
lambda = 2; mu = 1; Xzero = 1; % problem parameters
T = 1; N = 2^8; dt = 1/N;
dW = sqrt(dt)*randn(1,N); % Brownian increments
W = cumsum(dW); % discretized Brownian path
Xtrue = Xzero*exp((lambda-0.5*mu^2)*([dt:dt:T])+mu*W);
plot([0:dt:T],[Xzero,Xtrue],'m-'), hold on
R = 4; Dt = R*dt; L = N/R; % L EM steps of size Dt = R*dt
Xem = zeros(1,L); % preallocate for efficiency
Xtemp = Xzero;
for j = 1:L
Winc = sum(dW(R*(j-1)+1:R*j));
Xtemp = Xtemp + Dt*lambda*Xtemp + mu*Xtemp*Winc;
Xem(j) = Xtemp;
end
plot([0:Dt:T],[Xzero,Xem],'r--*'), hold off
xlabel('t','FontSize',12)
ylabel('X','FontSize',16,'Rotation',0,'HorizontalAlignment','right')
emerr = abs(Xem(end)-Xtrue(end))
%% 
close all 
clear 
clc

% initial condition
r=0.033;  L=0.2;  T = 1;  D = 4;    v=1/T;  DT=D*T;
omega1=v/r;  omega2=v/r;
 
N = 2^11;	dt = T/N;     
R = [1; 16; 32; 64; 128];
for p = 1:5
Dt = R(p)*dt; P = N/R(p);% Large Q--fast speed low accuracy

% initiallize
X_str = zeros(N,P);         Y_str = zeros(N,P);   
theta_str=zeros(N,P);       
Winc1 = zeros(N,P);         Winc2 = zeros(N,P);  
% 2 Wiener processes
randn('state',400)

dW1=sqrt(dt)*randn(N,N);
dW2=sqrt(dt)*randn(N,N);

% Straight SDE(0,0)-(0,1)
% Banana distribution

% for j=1:N 
%     for i = 2:P
%         Winc1 = sum(dW1(j,Q*(i-1)+1:Q*i));
%         Winc2 = sum(dW2(j,Q*(i-1)+1:Q*i));
%         theta_str(j,i) = theta_str(j,i-1) + Dt*r*(omega1-omega2)/L + sqrt(D)*r*(Winc1-Winc2)/L;
%     end
% end


for j=1:N
    for i=2:P
        Winc1 = sum(dW1(j,R*(i-1)+1:R*i));
        Winc2 = sum(dW2(j,R*(i-1)+1:R*i));
        theta_str(j,i) = theta_str(j,i-1) + Dt*r*(omega1-omega2)/L + sqrt(D)*r*(Winc1-Winc2)/L;
        X_str(j,i) = X_str(j,i-1) + 0.5*r*(omega1+omega2)*cos(theta_str(j,i-1))*Dt + 0.5*r*sqrt(D)*cos(theta_str(j,i-1))*(Winc1+Winc2);
        Y_str(j,i) = Y_str(j,i-1) + 0.5*r*(omega1+omega2)*sin(theta_str(j,i-1))*Dt + 0.5*r*sqrt(D)*sin(theta_str(j,i-1))*(Winc1+Winc2);
    end
end

Xmil(:,p) = X_str(:,end); % store Milstein solution at t =1
Ymil(:,p) = Y_str(:,end);
end
Xref = Xmil(:,1); % Reference solution
Yref = Ymil(:,1);
Xerr=abs(Xmil(:,2:5) - repmat(Xref,1,4));
% err = sqrt((Xmil(:,2:5) - repmat(Xref,1,4)).^2+(Ymil(:,2:5)-repmat(Yref,1,4)).^2); % Error in each path
err = abs(Xmil(:,2:5) - repmat(Xref,1,4))+abs(Ymil(:,2:5)-repmat(Yref,1,4)); % Error in each path
Dtvals = dt*R(2:5); % Milstein timesteps used
figure;
loglog(Dtvals,mean(err),'b*-'), hold on
% loglog(Dtvals,Dtvals,'r--'), hold off % reference slope of 1
% axis([1e-2 1e-1 1e-4 1])
xlabel('\Delta t')
ylabel('Sample average of | X(T) - X_L |')
title('EMstrong','FontSize',10)
