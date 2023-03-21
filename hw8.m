clear;
clc;
%% ode45
c1=1;
c2=2;
u=1;
dxdt=@(t,x)[x(2);-x(1)-((-c2-c1)/2)*x(2)+0.5*x(3)+(c2/2)*x(4);x(4);x(1)+c2*x(2)-x(3)-c2*x(4)+u];
[t,X]=ode45(dxdt,[0 40],[0;0;0;0]);
figure;
subplot(2,1,1);plot(t,X(:,1));
xlabel('Time t');
ylabel('Solution q1');
title('Solution with ODE45');
subplot(2,1,2);plot(t,X(:,3));
xlabel('Time t');
ylabel('Solution q2');

%% check stablity
i=1;
x=zeros();
for c1=-10:0.5:10
    for c2=-10:0.5:10

    A=[0 1 0 0;
       -1 (-c2-c1)/2 0.5 c2/2;
       0 0 0 1;
       1 c2 -1 -c2];
    B=[0;0;0;1];
    C=[1 0 0 0;
       0 0 1 0];
    D=[0;0];
   
    g=ss(A,B,C,D);  
    [num,den]=ss2tf(A,B,C,D);
    s=roots(den);
    if (s(1)<0)&&(s(2)<0)&&(s(3)<0)&&(s(4)<0)
        x(i,1)=c1;
        x(i,2)=c2;
        i=i+1;
    end
    end
end
%% check observability
c1=1;
c2=2;
 A=[0 1 0 0;
       -1 (-c2-c1)/2 0.5 c2/2;
       0 0 0 1;
       1 c2 -1 -c2];
    B=[0;0;0;1];
    C=[1 0 0 0;
       0 0 1 0];
    D=[0;0];
O=obsv(A,C);
or=rank(O)

%% check controbility
i=1;
x=zeros();
for c1=-10:0.5:10
    for c2=-10:0.5:10

    A=[0 1 0 0;
       -1 (-c2-c1)/2 0.5 c2/2;
       0 0 0 1;
       1 c2 -1 -c2];
    B=[0;0;0;1];
    C=[1 0 0 0;
       0 0 1 0];
    D=[0;0];
   Q=ctrb(A,B);
    if rank(Q)==4
        x(i,1)=c1;
        x(i,2)=c2;
        i=i+1;
    end
    end
end
%% check the controllability when apply u on mass 1
c1=1;
c2=2;
 A=[0 1 0 0;
       -1 (-c2-c1)/2 0.5 c2/2;
       0 0 0 1;
       1 c2 -1 -c2];
    B=[0;0.5;0;0];
    C=[1 0 0 0;
       0 0 1 0];
    D=[0;0];
 Q=ctrb(A,B);
 oq=rank(Q)
%% impulse
c1=1;
c2=2;
 A=[0 1 0 0;
       -1 (-c2-c1)/2 0.5 c2/2;
       0 0 0 1;
       1 c2 -1 -c2];
    B=[0;0;0;1];
    C=[1 0 0 0;
       0 0 1 0];
    D=[0;0];
    g=ss(A,B,C,D);  
    [num,den]=ss2tf(A,B,C,D);
    impulse(g);