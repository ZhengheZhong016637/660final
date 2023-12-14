clear all;
itermax = 500;
iter = 1;
tol = 1e-9;
axis = 1:itermax;
g_norm = zeros(itermax,1);
steps = zeros(itermax,182);
c = 0.1;
rho = 0.9;
A = readtable("Adjacency_matrix.csv");
A = A{:,:};
%X = [x;y]
X = normrnd(0,91,182,1);
f = potential(X,A);
steps(1,:) = X;
%gradient of potential is -force
g = -forces(X,A);
g_norm(1) = log(norm(g));
disp(log(norm(g)));
B = eye(182);
m = 10;

while iter<itermax
    p = -B\g;
    a=1;
    f_temp = potential(X+a*p,A);
    cpg = c*p'*g;
    while f_temp > f + a*cpg % check Wolfe's condition 1
        a = a*rho;
        if a < 1e-14
            fprintf("line search failed\n");
            fail_flag = 1;
            break;
        end
        f_temp = potential(X + a*p,A);        
    end
    iter = iter+1;
    s = a*p;
    X = X+s;    
    gnew = -forces(X,A);
    g_norm(iter) = log(norm(gnew));
    disp(log(norm(gnew)));
    f = potential(X,A);        
    y = gnew-g;
    B = Bupdate(B,y,s,iter,m);    
    g = gnew;    
    steps(iter,:) = X;
end

plot_graph(X,A)
figure;
plot(g_norm)
error = zeros(itermax,1);
for i = 1:itermax
    error(i) = norm(steps(i,:)-X);
end
figure;
plot(error);

function plot_graph(X,A)
    x = X(1:91);
    y = X(92:182);
    figure;
    hold on
    plot(x,y,'o','Markersize',15,'MarkerEdgeColor',[0.5,0,0],'MarkerFaceColor',[1,0,0]);
    ind = find(A == 1);
    [I,J] = ind2sub(size(A),ind);
    for k = 1 : length(ind)
        plot([x(I(k)),x(J(k))],[y(I(k)),y(J(k))],'linewidth',4,'Color',[0,0,0.5]);
    end
    daspect([1,1,1])
    axis off
end

function f = forces(X,A)
    x = X(1:91);
    y = X(92:182);
    % f = force = - grad U = column vector with 2*N components
    % x, y are column vectors with N components
    % A is an N-by-N adjacency matrix
    N = length(x);
    %% find pairwise distances between linked vertices
    xaux = x*ones(size(x))';
    yaux = y*ones(size(y))';
    dx = A.*xaux - A.*(xaux'); 
    dy = A.*yaux - A.*(yaux');
    dxy = sqrt(dx.^2+dy.^2);
    
    %% spring forces due to linked vertices
    Aind = find(A == 1);
    idiff = zeros(N);
    idiff(Aind) = 1 - 1./dxy(Aind);
    fx = -sum(idiff.*dx,2);
    afx = min(abs(fx),1);
    sfx = sign(fx);
    fx = afx.*sfx;
    
    fy = -sum(idiff.*dy,2);
    afy = min(abs(fy),1);
    sfy = sign(fy);
    fy = afy.*sfy;
    
    f_linked = [fx;fy];
    
    %% repelling spring forces due to unlinked vertices
    h = sqrt(3);
    Aind = find(A==0);
    A = ones(size(A))-A;
    dx = A.*xaux - A.*(xaux'); 
    dy = A.*yaux - A.*(yaux');
    dxy = sqrt(dx.^2+dy.^2);
    fac = zeros(N);
    diff = dxy - h;
    fac(Aind) = min(diff(Aind),0); 
    fx = sum(fac.*dx,2);
    fy = sum(fac.*dy,2);
    f_unlinked = -[fx;fy];
    
    f = f_linked + f_unlinked;
end

function U = potential(X,A)
    x = X(1:91);
    y = X(92:182);
    %linked potentials       
    xaux = x*ones(size(x))';
    yaux = y*ones(size(y))';
    dx = A.*xaux - A.*(xaux'); 
    dy = A.*yaux - A.*(yaux');
    dxy = 0.5*(sqrt(dx.^2+dy.^2)-1).^2;
    S = sum(dxy,"all");
    %unlinked potentials
    A = ones(size(A))-A;
    xaux = x*ones(size(x))';
    yaux = y*ones(size(y))';
    dx = A.*xaux - A.*(xaux'); 
    dy = A.*yaux - A.*(yaux');
    dxy = 0.5*(min(sqrt(dx.^2+dy.^2)-1,0)).^2;
    U = S+sum(dxy,"all");
end

function Bnew = Bupdate(B,y,s,iter,m)
    if mod(iter,m) == 0
        Bnew = eye(size(B));
    else
        Bnew = B-((B*s)*(s'*B))/(s'*B*s)+(y*y')/(y'*s);
    end
end