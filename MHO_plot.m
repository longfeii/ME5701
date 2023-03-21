
clear
clc
% initial parameters
l=0.2; r=0.033; v=1; %set the velocity
w1=v/r; w2=v/r;%two omega are same
ti=10000;% times(how many pathes)
T=1; %total time
dt=0.001; %time step
D=7; %noise coeffcient
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
        x(i,j) = x(i,j-1) + 0.5*r*(w1+w2)*cos(theta(i,j-1))*dt + sqrt(D)*0.5*r*cos(theta(i,j-1))*(dw1(i,j-1)+dw2(i,j-1))-0.5*sqrt(D)*0.5*r*cos(theta(i,j-1))*sqrt(D)*0.5*r*sin(theta(i,j-1))*(dw1(i,j-1)^2+dw2(i,j-1)^2-dt);
        y(i,j) = y(i,j-1) + 0.5*r*(w1+w2)*sin(theta(i,j-1))*dt + sqrt(D)*0.5*r*sin(theta(i,j-1))*(dw1(i,j-1)+dw2(i,j-1))+0.5*sqrt(D)*0.5*r*sin(theta(i,j-1))*sqrt(D)*0.5*r*cos(theta(i,j-1))*(dw1(i,j-1)^2+dw2(i,j-1)^2-dt);
        theta(i,j) = theta(i,j-1) + dt*r*(w1-w2)/l + sqrt(D)*r*(dw1(i,j-1)-dw2(i,j-1))/l;
    end
    plot(x(i,P),y(i,P),'b.')%plot every end point of evry path
    hold on
end

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
GMModel = fitgmdist(data,100);
GMMpdf=reshape(pdf(GMModel, [xc(:),yc(:)]),num,num);
hold on
contour(xc, yc, GMMpdf,3,'g','linewidth',2);
hold off
%% arc
clear
clc
% initial parameters
l=0.2; r=0.033; v=1; %set the velocity
a=1;ad=pi/2;
w1=ad*(a+l/2)/r; w2=ad*(a-l/2)/r;%two omega are different
ti=10000;% times(how many pathes)
T=1; %total time
dt=0.001; %time step
D=4; %noise coeffcient
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
        x(i,j) = x(i,j-1) + 0.5*r*(w1+w2)*cos(theta(i,j-1))*dt + sqrt(D)*0.5*r*cos(theta(i,j-1))*(dw1(i,j-1)+dw2(i,j-1))-0.5*sqrt(D)*0.5*r*cos(theta(i,j-1))*sqrt(D)*0.5*r*sin(theta(i,j-1))*(dw1(i,j-1)^2+dw2(i,j-1)^2-dt);
        y(i,j) = y(i,j-1) + 0.5*r*(w1+w2)*sin(theta(i,j-1))*dt + sqrt(D)*0.5*r*sin(theta(i,j-1))*(dw1(i,j-1)+dw2(i,j-1))+0.5*sqrt(D)*0.5*r*sin(theta(i,j-1))*sqrt(D)*0.5*r*cos(theta(i,j-1))*(dw1(i,j-1)^2+dw2(i,j-1)^2-dt);
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

