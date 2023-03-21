clear
clc
% initial parameters
l=0.2; r=0.033; v=1; %set the velocity
a=1;ad=pi/2;
w1=ad*(a+l/2)/r; w2=ad*(a-l/2)/r;%two omega are different
ti=10000;% times(how many pathes)
T=1; %total time
dt=0.001; %time step
D=1; %noise coeffcient
P=T/dt; %how many points in one path
x = zeros(ti,P);      
y = zeros(ti,P);   
theta=zeros(ti,P);

%Brownian increments
% randn('state',400)
dw1=sqrt(dt)*randn(ti,P);% two omega increase randomly and differently
dw2=sqrt(dt)*randn(ti,P);

%plot the distribution in cartesian
figure
for i=1:ti %for the first path
    for j=2:P % from start time to end
        x(i,j) = x(i,j-1) + 0.5*r*(w1+w2)*cos(theta(i,j-1))*dt + sqrt(D)*0.5*r*cos(theta(i,j-1))*(dw1(i,j-1)+dw2(i,j-1));
        y(i,j) = y(i,j-1) + 0.5*r*(w1+w2)*sin(theta(i,j-1))*dt + sqrt(D)*0.5*r*sin(theta(i,j-1))*(dw1(i,j-1)+dw2(i,j-1));
        theta(i,j) = theta(i,j-1) + dt*r*(w1-w2)/l + sqrt(D)*r*(dw1(i,j-1)-dw2(i,j-1))/l;
    end
    plot(x(i,P),y(i,P),'b.')%plot every end point of evry path
    hold on
end

%ideal path
t_ideal = linspace(0,pi/2,ti);
x_ideal=a*(sin(t_ideal));
y_ideal=a-a*cos(t_ideal);
plot(x_ideal,y_ideal,'--k')   
plot(a,a,'r*') 
axis([-0.5 2 -0.5 1.5])
xlabel('X position')
ylabel('Y position')
title(['DT=',num2str(D)])
hold on

%Cartesian pdf
%mean
mean_ca=zeros(1,2);
mean_ca(1)=sum(x(:,end))/ti;%sum all the end point's x
mean_ca(2)=sum(y(:,end))/ti;
%covariance
cov_ca=zeros(2,2);
cov_ca(1,1)=sum((x(:,end)-mean_ca(1)).^2)/ti;
cov_ca(1,2)=sum((x(:,end)-mean_ca(1)).*(y(:,end)-mean_ca(2)))/ti;
cov_ca(2,1)=cov_ca(1,2);
cov_ca(2,2)=sum((y(:,end)-mean_ca(2)).^2)/ti;
num=100;
xlim = get(gca,'XLim');
ylim = get(gca,'YLim');
[xc,yc] = meshgrid(linspace(xlim(1),xlim(2),num)',linspace(ylim(1),ylim(2),num)');
% zcc = reshape(mvnpdf([xc(:),yc(:)], mean_ca, [cov_ca(1,1),cov_ca(2,2)]),num,num);%get the two degree guassian value for each point
zc=zeros(num);
for i=1:num
    for j= 1:num
         zc(i,j)=exp((-0.5)*([xc(i,j),yc(i,j)]-mean_ca)*inv(cov_ca)*([xc(i,j),yc(i,j)]'-mean_ca'))/(2*pi*(det(cov_ca))^(0.5));
    end
end
contour(xc,yc,zc,3,'r','linewidth',2);
% %Exponential pdf
xx=x(:,end);     yy=y(:,end);
data=[xx, yy];
GMModel = fitgmdist(data,10);
GMMpdf=reshape(pdf(GMModel, [xc(:),yc(:)]),num,num);
hold on
c2=contour(xc, yc, GMMpdf,2,'g','linewidth',2);
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
axis([1 2 -0.5 0.5])
xlabel('X position')
ylabel('Y position')
title('exponential coordinates')
%pdf in exponential coordinate
%mean
mean_exp=zeros(1,2);
mean_exp(1)=sum(x_exp)/ti;%sum all the end point's x
mean_exp(2)=sum(y_exp)/ti;
%covariance
cov_exp=zeros(2,2);
cov_exp(1,1)=sum((x_exp-mean_exp(1)).^2)/ti;
cov_exp(1,2)=sum((x_exp-mean_exp(1)).*(y_exp-mean_exp(2)))/ti;
cov_exp(2,1)=cov_exp(1,2);
cov_exp(2,2)=sum((y_exp-mean_exp(2)).^2)/ti;
xlim_exp = get(gca,'XLim');
ylim_exp = get(gca,'YLim');
[xe,ye] = meshgrid(linspace(xlim_exp(1),xlim_exp(2),num)',linspace(ylim_exp(1),ylim_exp(2),num)');
zc_exp=zeros(num);
for i=1:num
    for j= 1:num
         zc_exp(i,j)=exp((-0.5)*([xe(i,j),ye(i,j)]-mean_exp)*inv(cov_exp)*([xe(i,j),ye(i,j)]'-mean_exp'))/(2*pi*(det(cov_exp))^(0.5));
    end
end
contour(xe,ye,zc_exp,3,'r','linewidth',2);