function f = randomfunction(m)

gammas=normrnd(0,1,[1,m]);
c=(1:m);
if m>10
    c(1:10)=1;
    c(11:m)=ones(1,m-10)./(c(11:m)-10).^2;
else
    c(1:m)=1;
end
c=c.*gammas;
onedim = @(x) sum(c.*[cos(x.*(1:m/2)) sin(x.*(1:m/2))]);
f = @(x) arrayfun(onedim,x);
end

