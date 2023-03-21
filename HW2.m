clear
%a
X=[0;1;-1;-2];
Y=[0;1;0.5;4];
n=3;
A=fitbypolynomial(X,Y,n);
At=zeros(size(A));
for i=1:n+1
At(i)=A(n+2-i);
end
display(At)
x=-2:0.01:1;
y=polyval(At,x);
plot(X,Y,'o');
hold on;
plot(x,y,'r');
%b
A2=fitbypseudo(X,Y,n);
A2t=zeros(size(A2));
for i=1:n
A2t(i)=A2(n+1-i);     
end
display(A2t)
y2=polyval(A2t,x);
hold on;
plot(x,y2,'g');
%c
syms x a
y=int((a.*x.^4+(2*a-0.25).*x.^3+(0.75-a).*x.^2+(0.5-2*a).*x).^2,x,-2,1);
a=(1791/560)/(2*81/70);
A4=[a,2*a-0.25,0.75-a,0.5-2*a,0]
x=-2:0.01:1;
y3=polyval(A4,x);
hold on;
plot(x,y3,'b');
legend('dot','a','b','c');
%% 
x=-2:0.01:1;
A5=[a,b,c,d,e];
y4=
C=integral(y4^2,-2,1)



%% 
Atoda=diag(repmat([2],1,100))+diag(repmat([-1],1,99),1)+diag(repmat([-1],1,99),-1);
u=zeros(100,100);
u(:,1)=rand(100,1);
for k=1:100
    u(:,k+1)=Atoda*u(:,k)/(norm(Atoda*u(:,k)));
    c(k)=(u(:,k)')*(Atoda*u(:,k));
    e(k)=norm(c(k)*u(:,k)-Atoda*u(:,k));
end
figure;
plot(c);
figure;
plot(e);
%% 

%8
f1=fs(1,1);f2=fs(1,2);f3=fs(1,5);f4=fs(1,10);
figure;
subplot(2,2,1);fplot(f1);
subplot(2,2,2);fplot(f2);
subplot(2,2,3);fplot(f3);
subplot(2,2,4);fplot(f4);
function f=fs(h,B)
syms x
f=h/2;
for i=1:2:B
    f=f+((-4)*h/((i^2)*(pi^2)))*cos(i*x);
end
end

function [A]=fitbypolynomial(X,Y,n)
m=size(X);
A=zeros(m(1));
X2=zeros(m(1),m(1));
for i=1:n+1
for j=1:n+1
    X2(i,j)=(X(i))^(j-1);
end
end
    A=X2\Y;
   
end

function [A]=fitbypseudo(X,Y,n)
m=size(X);
A=zeros(m(1)-1);
X2=zeros(m(1),m(1)-1);
for i=1:n+1
for j=1:n
    X2(i,j)=(X(i))^(j-1);
end
end
    A=(inv(transpose(X2)*X2))*transpose(X2)*Y;
end


