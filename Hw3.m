clear
n=10;
x=linspace(0,10,n);
y=linspace(0,10,n);
z=rand(n,n);
[xi,yi]=meshgrid(x,y);
subplot(2,2,1);Zgrid=griddata(x,y,z,xi,yi,'v4');surf(xi,yi,Zgrid);title('fit all the point with grid data method');xlabel('X'), ylabel('Y'), zlabel('Z')

for i=1:10 A=fitbypolynomial(X,Y,n);
At=zeros(size(A));
for i=1:n+1
At(i)=A(n+2-i);
end
x=-2:0.01:1;
y=polyval(At,x);
subplot(2,2,2);
plot(x,y,'r');
hold on;title('fit all the point with pseudo-inverse method');xlabel('X'), ylabel('Y'), zlabel('Z')
% subplot(2,2,3);fplot(f3);title('fit all the point without method');xlabel('X'), ylabel('Y'), zlabel('Z')
% subplot(2,2,4);fplot(f4);
