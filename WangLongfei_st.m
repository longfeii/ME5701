clear
clc
% initial parameters
l=0.2; r=0.033; v=1; %set the velocity
w1=v/r; w2=v/r;%two omega are same
ti=10000;% times(how many pathes)
T=1; %total time
dt=0.001; %time step
D=4; %noise coeffcient
P=T/dt; %how many points in one path
x = zeros(ti,P);      
y = zeros(ti,P);   
theta=zeros(ti,P);

%Brownian increments
randn('state',400)
dw1=sqrt(dt)*randn(ti,P);% two omega increase randomly and differently
dw2=sqrt(dt)*randn(ti,P);

%plot the distribution in cartesian
figure %EM method of SDE 
for i=1:ti %for the first path
    for j=2:P % from start time to end
        x(i,j) = x(i,j-1) + 0.5*r*(w1+w2)*cos(theta(i,j-1))*dt + sqrt(D)*0.5*r*cos(theta(i,j-1))*(dw1(i,j-1)+dw2(i,j-1));
        y(i,j) = y(i,j-1) + 0.5*r*(w1+w2)*sin(theta(i,j-1))*dt + sqrt(D)*0.5*r*sin(theta(i,j-1))*(dw1(i,j-1)+dw2(i,j-1));
        theta(i,j) = theta(i,j-1) + dt*r*(w1-w2)/l + sqrt(D)*r*(dw1(i,j-1)-dw2(i,j-1))/l;
    end
    plot(x(i,P),y(i,P),'b.')%plot every end point of evry path
    hold on
end
save('Data',"x","y")
%ideal path
t_ideal = linspace(0,1,ti);
x_ideal=v*t_ideal;
y_ideal=zeros(1,ti);
plot(x_ideal,y_ideal,'--k')   
plot(1,0,'r*') 
axis([-0.5 1.5 -1 1])
xlabel('X position')
ylabel('Y position')
title(['DT=',num2str(D)])
hold on
%Cartesian pdf
%mean
mean_ca=zeros(1,3);
mean_ca(1)=sum(x(:,end))/ti;%sum all the end point's x
mean_ca(2)=sum(y(:,end))/ti;
mean_ca(3)=sum(theta(:,end))/ti;
%covariance
multi=zeros(3);
for o=1:ti
    multi=multi+([x(o,end)-mean_ca(1);y(o,end)-mean_ca(2);theta(o,end)-mean_ca(3)]*[x(o,end)-mean_ca(1);y(o,end)-mean_ca(2);theta(o,end)-mean_ca(3)]');
end
cov_ca = multi/ti;
num=100;
xlim = get(gca,'XLim');
ylim = get(gca,'YLim');
[xc,yc] = meshgrid(linspace(xlim(1),xlim(2),num)',linspace(ylim(1),ylim(2),num)');
zc=zeros(num);
for i=1:num
    for j= 1:num
         zc(i,j)=exp((-0.5)*([xc(i,j),yc(i,j)]-mean_ca(1:2))*inv(cov_ca(1:2,1:2))*([xc(i,j),yc(i,j)]'-mean_ca(1:2  )'))/(2*pi*(det(cov_ca(1:2,1:2)))^(0.5));
    end
end
c1=contour(xc,yc,zc,3,'r','linewidth',2);
% %Exponential pdf
xx=x(:,end);     yy=y(:,end);
data=[xx, yy];
GMModel = fitgmdist(data,10);
GMMpdf=reshape(pdf(GMModel, [xc(:),yc(:)]),num,num);
hold on
c2=contour(xc, yc, GMMpdf,2,'g','linewidth',2);
hold off
legend([c1,c2],{'Cart pdf','Exp pdf'})
%QQ plot for Cartesian
dm=zeros(1,ti);
for i=1:ti
    dm(i)=[x(i,end)-mean_ca(1),y(i,end)-mean_ca(2)]*inv(cov_ca(1:2,1:2))*[x(i,end)-mean_ca(1),y(i,end)-mean_ca(2)]';
end
dm=sort(dm);
xt=0:40;
yt=xt;
pt=((1:ti)-0.5)/ti;  
xm=chi2inv(pt,3); 
figure;
scatter(xm,dm','.');
hold on
plot(xt,yt);
hold off
%% Express points by exponential coordinate
x_exp=zeros(1,ti);
y_exp=zeros(1,ti);
a_exp=zeros(1,ti);
%plot distribution in exponential
figure;
for i=1:ti
    H=[cos(theta(i,end)),-sin(theta(i,end)),x(i,end);
       sin(theta(i,end)),cos(theta(i,end)),y(i,end);
        0,0,1 ];
    N=logm(H); %exp(N)=H 
    x_exp(i)=N(1,3);        % v1 
    y_exp(i)=N(2,3);        % v2
    a_exp(i)=N(2,1);    % alpha
    plot(x_exp(i),y_exp(i),'b.');
    hold on
end
axis([0.4 1.5 -0.8 0.8])
xlabel('X position')
ylabel('Y position')
title('exponential coordinates')
%pdf in exponential coordinate
%mean
mean_exp=zeros(1,3);
mean_exp(1)=sum(x_exp)/ti;%sum all the end point's x
mean_exp(2)=sum(y_exp)/ti;
mean_exp(3)=sum(a_exp(:,end))/ti;
%covariance
multi=zeros(3);
for o=1:ti
    multi=multi+([x_exp(o)-mean_exp(1);y_exp(o)-mean_exp(2);a_exp(o)-mean_exp(3)]*[x_exp(o)-mean_exp(1);y_exp(o)-mean_exp(2);a_exp(o)-mean_exp(3)]');
end
cov_exp = multi/ti;
xlim_exp = get(gca,'XLim');
ylim_exp = get(gca,'YLim');
[xe,ye] = meshgrid(linspace(xlim_exp(1),xlim_exp(2),num)',linspace(ylim_exp(1),ylim_exp(2),num)');
zc_exp=zeros(num);
for i=1:num
    for j= 1:num
         zc_exp(i,j)=exp((-0.5)*([xe(i,j),ye(i,j)]-mean_exp(1:2))*inv(cov_exp(1:2,1:2))*([xe(i,j),ye(i,j)]'-mean_exp(1:2)'))/(2*pi*(det(cov_exp(1:2,1:2)))^(0.5));
    end
end
contour(xe,ye,zc_exp,3,'r','linewidth',2); 
%QQ plot for exponential
dme=zeros(1,ti);
for i=1:ti
    dme(i)=[x_exp(i)-mean_exp(1),y_exp(i)-mean_exp(2)]*inv(cov_exp(1:2,1:2))*[x_exp(i)-mean_exp(1),y_exp(i)-mean_exp(2)]';
end
dme=sort(dme);
xt=0:40;
yt=xt;
pt=((1:ti)-0.5)/ti;  
xm=chi2inv(pt,3); 
figure;
scatter(xm,dme','.');
hold on
plot(xt,yt);
hold off
%% Propagation method
% mean
t=T;
mean_prop=[1 0 r*w1*t; 0 1 0; 0 0 1];
% covariance
sigma11=0.5*D*t*r^2;
sigma22=(2*D*(w1^2)*(r^4)*(t^3))/( 3*(l^2) );     
sigma23=D*w1*r^3*t^2/l^2;
sigma32=sigma23;
sigma33=2*D*r^2*t/(l^2);
cov_prop=[sigma11 0 0; 0 sigma22 sigma23; 0 sigma32 sigma33];
% zc_exp2=zeros(num);
% for i=1:num
%     for j= 1:num
%          zc_exp2(i,j)=exp((-0.5)*([xe(i,j),ye(i,j)]-mean_prop(1:2))*inv(cov_prop(1:2,1:2))*([xe(i,j),ye(i,j)]'-mean_prop(1:2)'))/(2*pi*(det(cov_prop(1:2,1:2)))^(0.5));
%     end
% end
% contour(xe,ye,zc_exp2,3,'g','linewidth',2); 
%% Kalman filter
kf = [];
xk=[0;0];
Pk=zeros(2);
A=[1,0;0,0];
R=[1,0;
   0,1];
H=eye(2);
% Q=[0.1,0;
%    0,0.1];
for i = 1:P
    z= [x(1,i);y(1,i)];
    N1=[sqrt(D)*0.5*r*cos(theta(1,i)),sqrt(D)*0.5*r*cos(theta(1,i));
       sqrt(D)*0.5*r*sin(theta(1,i)),sqrt(D)*0.5*r*sin(theta(1,i))];
    Q=N1*N1';
    % prediction
    xhat = A * xk;
    Phat = A * Pk * A' + Q;
    %gain
    K = Phat * H' * inv( H*Phat*H'+R);
    xhat = xhat + K*(z-H*xhat);
    Pk = (1 - K * H) * Phat;
    % Store
    xk=xhat;
    kf = [kf,xhat];
end
figure
plot(x(1,:),y(1,:),'linewidth', 1,'color','r');
hold on
plot(kf(1,:),kf(2,:),'b');
plot(x_ideal,y_ideal,'--k');
legend('Actual path','KF estimation path','ideal path');
hold off;