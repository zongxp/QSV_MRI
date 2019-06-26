function x=randn_white(n)
% x=randn_white(n)
tic;
if length(n)==1
    n=[1,n];
end
x=zeros(n);
for i=1:length(x(:))
s =2;
while s>=1
U1 = rand;
U2 = rand;
V1=2*U1-1;
V2 = 2*U2-1;
s = V1*V1+V2*V2;
end

x(i) = sqrt(-2*log(s)/s)*V1;
% if mod(i,10000)==0
%     
%     disp(i);
%     disp(toc);
%     tic;
% end
end



   