%% 3.random walk - Irene
N = 100;% Amount of random walks
M = 10; % The length of the random walk
g=(1/4)^M; % the probability of random walk 
z=zeros(1,N);
w=zeros(1,N);
for n=1:N
    x_t = zeros(1,M);
    y_t = zeros(1,M);
    x_t(1) = 0;
    y_t(1) = 0;
    for m=1:(M-1)
        A=randi(4);
        if A==1
            x_t(m+1)= x_t(m) + 1;
            y_t(m+1)= y_t(m);
        elseif A==2
            x_t(m+1)= x_t(m) - 1;
            y_t(m+1)= y_t(m);
        elseif A==3
            x_t(m+1)= x_t(m) ;
            y_t(m+1)= y_t(m) + 1 ;
        else
            x_t(m+1)= x_t(m) ;
            y_t(m+1)= y_t(m) - 1 ;
        end
    end
    x_t;
    y_t;   
    % When is it a random walk?
    z_t=1;
    for i=1:M 
        for j=setdiff(1:M, i) 
            if ((y_t(i) == y_t(j)) && (x_t(j)==x_t(i)))
                z_t=0;
            end
        end
    end
    z(n)=z_t;
    w_t= (z_t)/g;
    w(n)=w_t;
end     
z;
s=mean(z);
ww=s/g;
N_SA=sum(z);
c_n = (1/N)*sum(ww)
c_na = (1/N)*N_SA*(4^M) % right way
%% 

%% 4.Sequential improtance sampling. Notations here will try to follow PPT notations.

N = 5000;% Number of simulated particles
% part = zeros(N,2); % init all particles, started at (0,0) in Z2 dimension. Column 1 represents for the x axis, and column 2 for y axis.
n = 200;
w = zeros(n+1,N); % init w. 
w(1,:) = ones(1,N);
RW = zeros(n+1,2,N); % Use multidimensional array to store generated Random Walks. And init all particles, started at (0,0) in Z2 dimension. Column 1 represents for the x axis, and column 2 for y axis, pages for different particles.
cn = zeros(n,1); % Init cn. It is the value we want. Alternative: usedouble.empty(0,n)

for k = 1:n
    for i = 1:N
        next = getNext(RW(1:k,:,i));
        possibles = size(next,1);
        if possibles == 0
            gk1 = 0; %Here gk1 means gk+1 in slides
            wk1 = 0;
        else
            gk1 = 1.0 / possibles;
            wk1 = (1.0 / gk1) * w(k,i);
            U = rand(1);
            xk1 = next(floor(1+possibles*U),:);
            RW(k+1,:,i) = xk1;
        end
        w(k+1,i) = wk1;
    end
    cn(k) = mean(w(k+1,:));
end

plot(cn)% estimation of cn growing with n
j = 1 % choose the jth particle 
plot(RW(:,1,j),RW(:,2,j)) % plot the randon walk trace
%% 5. SISR
N = 5000;% Number of simulated particles
% part = zeros(N,2); % init all particles, started at (0,0) in Z2 dimension. Column 1 represents for the x axis, and column 2 for y axis.
n = 200;
w = zeros(n+1,N); % init w. 
w(1,:) = ones(1,N);
RW = zeros(n+1,2,N); % Use multidimensional array to store generated Random Walks. And init all particles, started at (0,0) in Z2 dimension. Column 1 represents for the x axis, and column 2 for y axis, pages for different particles.
cn = zeros(n,1); % Init cn. It is the value we want. Alternative: usedouble.empty(0,n)

for k = 1:n
    for i = 1:N
        next = getNext(RW(1:k,:,i));
        possibles = size(next,1);
        if possibles == 0
            gk1 = 0; %Here gk1 means gk+1 in slides
            wk1 = 0;
        else
            gk1 = 1.0 / possibles;
            wk1 = (1.0 / gk1) * w(k,i);
            U = rand(1);
            xk1 = next(floor(1+possibles*U),:);
            RW(k+1,:,i) = xk1;
        end
        w(k+1,i) = wk1;
    end
    cn(k) = mean(w(k+1,:));
end