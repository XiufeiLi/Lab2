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
j = 100 % choose the jth particle 
plot(RW(:,1,j),RW(:,2,j)) % plot the randon walk trace
% Since I restore whole w, it's possible to plot w distribution in each
% step with function histogram
% Below is for problem 5 plot, must use w in problem 4, since w will be
% changed in problem 5.
k = 150; % choose the k-th step u want to see details
CW = cumsum([0 w(k,:)]);
[num,ind] = histc(rand(1,N), CW/CW(end));
plot(num) % after resampling, the replication numbers of each particle.
plot(sort(ind)) % For debug
%% 5. SISR
N = 5000;% Number of simulated particles
% part = zeros(N,2); % init all particles, started at (0,0) in Z2 dimension. Column 1 represents for the x axis, and column 2 for y axis.
n = 200;
w = zeros(n+1,N); % init w. 
w(1,:) = ones(1,N);
RW = zeros(n+1,2,N); % Use multidimensional array to store generated Random Walks. And init all particles, started at (0,0) in Z2 dimension. Column 1 represents for the x axis, and column 2 for y axis, pages for different particles.
cn = zeros(n+1,1); % Init cn. It is the value we want. Alternative: usedouble.empty(0,n) Here n+1 means that the first one is init state. Cause here we need recursive calculate the cn based on previous cn, so add a init state for simplicity.
cn(1) = 1;
for k = 1:n
    for i = 1:N
        next = getNext(RW(1:k,:,i));
        possibles = size(next,1);
        if possibles == 0
            gk1 = 0; %Here gk1 means gk+1 in slides
            wk1 = 0;
        else
            gk1 = 1.0 / possibles;
            wk1 = (1.0 / gk1);
            U = rand(1);
            xk1 = next(floor(1+possibles*U),:);
            RW(k+1,:,i) = xk1;
        end
        w(k+1,i) = wk1;
    end
%     mean(w(k+1,:))
    cn(k+1) = cn(k)*mean(w(k+1,:));
    % selection part
    CW = cumsum([0 w(k+1,:)]);
    [~,ind] = histc(rand(1,N), CW/CW(end));
    RW(:,:,:) = RW(:,:,ind);
end
% Using result of 5 for problem 6. Only one time calculation for test.
x = 1:n;
X = [x',log(x)'];
y = log(cn(2:end));
md = fitlm(X, y, 'Intercept', true); %Reference: https://www.mathworks.com/help/stats/linear-regression-model-workflow.html
coef = md.Coefficients.Estimate;
A2 = exp(coef(1))
mu2 = exp(coef(2))
gamma2 = coef(3) + 1
%% Problem 6. It needs several iteration of problem 5 and runs slow, so I stored the value in annotation at bottom.
T = 20; % Estimation times
A2 = zeros(T,1);
mu2 = zeros(T,1);
gamma2 = zeros(T,1);
for t = 1:T % 20 times estimation
    N = 5000;% Number of simulated particles
    % part = zeros(N,2); % init all particles, started at (0,0) in Z2 dimension. Column 1 represents for the x axis, and column 2 for y axis.
    n = 200;
    w = zeros(n+1,N); % init w. 
    w(1,:) = ones(1,N);
    RW = zeros(n+1,2,N); % Use multidimensional array to store generated Random Walks. And init all particles, started at (0,0) in Z2 dimension. Column 1 represents for the x axis, and column 2 for y axis, pages for different particles.
    cn = zeros(n+1,1); % Init cn. It is the value we want. Alternative: usedouble.empty(0,n) Here n+1 means that the first one is init state. Cause here we need recursive calculate the cn based on previous cn, so add a init state for simplicity.
    cn(1) = 1;
    for k = 1:n
        for i = 1:N
            next = getNext(RW(1:k,:,i));
            possibles = size(next,1);
            if possibles == 0
                gk1 = 0; %Here gk1 means gk+1 in slides
                wk1 = 0;
            else
                gk1 = 1.0 / possibles;
                wk1 = (1.0 / gk1);
                U = rand(1);
                xk1 = next(floor(1+possibles*U),:);
                RW(k+1,:,i) = xk1;
            end
            w(k+1,i) = wk1;
        end
%         mean(w(k+1,:))
        cn(k+1) = cn(k)*mean(w(k+1,:));
        % selection part
        CW = cumsum([0 w(k+1,:)]);
        [~,ind] = histc(rand(1,N), CW/CW(end));
        RW(:,:,:) = RW(:,:,ind);
    end
    % Using result of 5 for problem 6. Only one time calculation for test.
    x = 1:n;
    X = [x',log(x)']; % ln transfer, described in problem statement.
    y = log(cn(2:end));
    md = fitlm(X, y, 'Intercept', true); %Reference: https://www.mathworks.com/help/stats/linear-regression-model-workflow.html
    coef = md.Coefficients.Estimate;
    A2(t) = exp(coef(1)); % Transfer back.
    mu2(t) = exp(coef(2));
    gamma2(t) = coef(3) + 1;
end
% Results:
% A2 = [1.45783776187533,1.49661616643717,1.34579370173985,1.35408772286713,1.37054844335199,1.34953753399444,1.34336240362084,1.34175139324528,1.32659437605983,1.41201740735313,1.42243400326149,1.33565332859375,1.29019975639461,1.29516889274759,1.51446003601618,1.42341345347709,1.32012078198599,1.35495736766942,1.42149197086665,1.41000489774257]
% mu2 = [2.64033365862045,2.64268398629262,2.64102730283201,2.63917147349139,2.63667118296124,2.63736210872200,2.63581480044301,2.63655908612701,2.64251571650845,2.64423781810494,2.64058018933594,2.63563417373426,2.63882696831863,2.64036088640296,2.64126254064595,2.63909919632319,2.63824534266435,2.64089706890159,2.63553199765645,2.64045679321556]
% gamma2 = [1.27380717439400,1.25275243032653,1.30208769222917,1.30478377683079,1.30574303971811,1.31974785654294,1.32320024842473,1.32273262168226,1.31418397846196,1.27791810298498,1.28504285819992,1.32132029987453,1.33855649378919,1.31903721157300,1.25451469019331,1.29167294497936,1.32541321349035,1.29859792924225,1.29667143606282,1.28974922008958]
%% problem 9. d dimension. SISR as problem 5, coefficients as problem 6. Similarly, coefficients estimation results are in bottom annotation.
d = 3; % Dimension
T = 20; % Estimation times
Ad = zeros(T,1);
mud = zeros(T,1);
gammad = zeros(T,1);
for t = 1:T % 20 times estimation
    N = 5000;% Number of simulated particles
    % part = zeros(N,2); % init all particles, started at (0,0) in Z2 dimension. Column 1 represents for the x axis, and column 2 for y axis.
    n = 200;
    w = zeros(n+1,N); % init w. 
    w(1,:) = ones(1,N);
    RW = zeros(n+1,d,N); % Use multidimensional array to store generated Random Walks. And init all particles, started at (0,0) in Z2 dimension. Column 1 represents for the x axis, and column 2 for y axis, pages for different particles.
    cn = zeros(n+1,1); % Init cn. It is the value we want. Alternative: usedouble.empty(0,n) Here n+1 means that the first one is init state. Cause here we need recursive calculate the cn based on previous cn, so add a init state for simplicity.
    cn(1) = 1;
    for k = 1:n
        for i = 1:N
            next = getNextd(RW(1:k,:,i));
            possibles = size(next,1);
            if possibles == 0
                gk1 = 0; %Here gk1 means gk+1 in slides
                wk1 = 0;
            else
                gk1 = 1.0 / possibles;
                wk1 = (1.0 / gk1);
                U = rand(1);
                xk1 = next(floor(1+possibles*U),:);
                RW(k+1,:,i) = xk1;
            end
            w(k+1,i) = wk1;
        end
%         mean(w(k+1,:))
        cn(k+1) = cn(k)*mean(w(k+1,:));
        % selection part
        CW = cumsum([0 w(k+1,:)]);
        [~,ind] = histc(rand(1,N), CW/CW(end));
        RW(:,:,:) = RW(:,:,ind);
    end
    % Using result of 5 for problem 6. Only one time calculation for test.
    x = 1:n;
    X = [x',log(x)'];
    y = log(cn(2:end));
    md = fitlm(X, y, 'Intercept', true); %Reference: https://www.mathworks.com/help/stats/linear-regression-model-workflow.html
    coef = md.Coefficients.Estimate;
    Ad(t) = exp(coef(1));
    mud(t) = exp(coef(2));
    gammad(t) = coef(3) + 1;
end
mudBound = 2*d - 1 - 1/(2d) - 3/(2d)^2 - 16/(2d)^3 % Asymptotic bound on mud for large d found in Graham (2014);