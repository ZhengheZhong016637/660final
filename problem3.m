clear all;
beta = 0.2:0.01:1;
mean = zeros(81,1);
var = zeros(81,1);
miu = zeros(81,1);
for k = 1:81    
    if beta(k)<0.4408
        miu(k) = 0;
    else
        miu(k) = (1-(sinh(2*beta(k)))^(-4))^(1/8);
    end
    [mean(k),var(k)] = mu(beta(k));    
end
var1 = sqrt(var)+mean;
var2 = -sqrt(var)+mean;
figure;
hold on
a1 = plot(1:81,mean); M1 = "numerical mean";
a2 = plot(1:81,miu); M2 = "analytic value";
a3 = plot(1:81,var1,"--"); M3 = "mean + sqrt(var)";
a4 = plot(1:81,var2,"--"); M4 = "mean - sqrt(var)";
legend([a1,a2,a3,a4], [M1, M2,M3,M4]);

function [mean,var] = mu(beta)
    s = ones(30,30);
    maxiter = 10^8;
    m = sum(s,"all")/900;
    mean = m;
    var = 0;
    for iter = 1:maxiter        
        i = randi([1 30]);
        j = randi([1 30]);
        deltah = deltaH(s,i,j);
        if deltah<0
            accept = true;
        else
            u = rand(1);
            if u<exp(-beta*deltah)
                accept = true;
            else
                accept = false;
            end
        end
        if accept 
            s(i,j) = -s(i,j);
            m = sum(s,"all")/900;
        end
        mean = (iter*mean+m)/(iter+1);
        var = ((iter-1)*var+(m-mean)^2)/iter;
    end
end

function m = H(s)
    m=0;
    for i = 1:30
        for j = 1:30
            if i+1 == 31
                k = 1;
            else 
                k = i+1;
            end
            if j+1 == 31
                l = 1;
            else
                l = j+1;
            end
            m = m-s(i,j)*(s(k,j)+s(i,l));
        end
    end
end

function d = deltaH(s,i,j)   
    ds = s(i,j)-(-1)*s(i,j);
    if i+1 == 31
        k = 1;
    else 
        k = i+1;
    end
    if j+1 == 31
        l = 1;
    else
        l = j+1;
    end
    if i-1 == 0
        a = 30;
    else
        a = i-1;
    end
    if j-1 == 0
        b = 30;
    else
        b = j-1;
    end
    d = ds*(s(k,j)+s(i,l)+s(a,j)+s(i,b));
end