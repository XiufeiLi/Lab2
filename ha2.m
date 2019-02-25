%% 3.random walk - Irene
N = 5000;% Amount of random walks
M = 200; % The length of the random walk
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
%% 3.

N = 5000;% Number of simulated particles
n = 200;
w = zeros(n+1,N); % init w. 
w(1,:) = ones(1,N);
RW = zeros(n+1,2,N); % Use multidimensional array to store generated Random Walks. And init all particles, started at (0,0) in Z2 dimension. Column 1 represents for the x axis, and column 2 for y axis, pages for different particles.
cn = zeros(n,1); % Init cn. It is the value we want. Alternative: usedouble.empty(0,n)

for k = 1:n
    for i = 1:N
        next = getNextSRW(RW(k,:,i)); % Possible next positions
        possibles = size(next,1);
        U = rand(1);
        xk1 = next(floor(1+possibles*U),:);
        RW(k+1,:,i) = xk1;
        gk1 = 1.0 / possibles;
        if w(k,i) == 0
            wk1 = 0;
        else
            whetherSAW(RW(1:k+1,:,i));
            if whetherSAW(RW(1:k+1,:,i))
                wk1 = (1.0 / gk1) * w(k,i);
            else
                wk1 = 0;
            end
        end
        wk1;
        RW(k+1,:,i) = xk1;
        w(k+1,i) = wk1;
    end
    cn(k) = mean(w(k+1,:));
end
plot(cn)
xlabel("n")
ylabel("C_n(2)")
title("Estimation C_n(2) with 5000 paticles for 200 steps, SIS")
set(gcf,'unit','centimeters','position',[10 5 20 12])
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
        next = getNext(RW(1:k,:,i)); % Possible next positions
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
plot(log(cn))
xlabel("n")
ylabel("ln(C_n(2))")
title("Estimation C_n(2) with 5000 paticles for 200 steps, SIS based on g_n")
set(gcf,'unit','centimeters','position',[10 5 20 12])

j = 150 % choose the jth particle 
index = find(RW(:,1,j) == 0 & RW(:,2,j) == 0);
endpoint = index(2)
plot(RW(1:endpoint-1,1,j),RW(1:endpoint-1,2,j)) % plot the randon walk trace
set(gcf,'unit','centimeters','position',[10 5 10 10])
xlabel("x")
ylabel("y")
title("Generated random walk")

bar(w(201,:),3)
set(gcf,'unit','centimeters','position',[10 5 20 10])
xlim([0,5000])
xlabel("Particle index")
ylabel("Weight value")
title("Weights at step 200")

% histogram(w(11,:),500)
% set(gcf,'unit','centimeters','position',[10 5 20 10])
% set(gca,'xscale','log')
% xlim([1e-15,1e5])


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
        next = getNext(RW(1:k,:,i)); % Possible next positions
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

plot(cn(2:end))
xlabel("n")
ylabel("C_n(2)")
xlim([0,210])
title("Estimation C_n(2) with 5000 paticles for 200 steps, SISR")
set(gcf,'unit','centimeters','position',[10 5 20 12])

plot(log(cn(2:end)))
xlabel("n")
ylabel("ln(C_n(2))")
xlim([0,200])
title("Estimation C_n(2) with 5000 paticles for 200 steps, SISR")
set(gcf,'unit','centimeters','position',[10 5 20 12])

j = 100 % choose the jth particle 
index = find(RW(:,1,j) == 0 & RW(:,2,j) == 0);
plot(RW(:,1,j),RW(:,2,j)) % plot the randon walk trace
set(gcf,'unit','centimeters','position',[10 5 10 10])
xlabel("x")
ylabel("y")
title("Generated self-avoiding random walk")

bar(w(100,:))
set(gcf,'unit','centimeters','position',[10 5 20 10])
xlim([0,5000])
xlabel("Particle index")
ylabel("Weight value")
title("Weights at step 100, SISR")

histogram(w(201,:),6)
set(gcf,'unit','centimeters','position',[10 5 20 10])
xlabel("Weights")
ylabel("Numbers")
title("Weights at step 200, SISR")
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
            next = getNext(RW(1:k,:,i)); % Possible next positions
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
plot(gamma2)
xlabel("Estimations index")
ylabel("\gamma_2")
title("Estimation of \gamma_2")
%% problem 9. d dimension. SISR as problem 5, coefficients as problem 6. Similarly, coefficients estimation results are in bottom annotation. It is super slow...
d = 3; % Dimension. Be careful with d = 4, which is not supported here, since it has another Cn(d) expression.
T = 20; % Estimation times
Ad = zeros(T,1);
mud = zeros(T,1);
gammad = zeros(T,1);
for t = 1:T % 20 times estimation
    N = 2000;% Number of simulated particles. Reduced to 2000 since out of memory.
    % part = zeros(N,2); % init all particles, started at (0,0) in Z2 dimension. Column 1 represents for the x axis, and column 2 for y axis.
    n = 200;
    w = zeros(n+1,N); % init w. 
    w(1,:) = ones(1,N);
    RW = zeros(n+1,d,N); % Use multidimensional array to store generated Random Walks. And init all particles, started at (0,0) in Z2 dimension. Column 1 represents for the x axis, and column 2 for y axis, pages for different particles.
    cn = zeros(n+1,1); % Init cn. It is the value we want. Alternative: usedouble.empty(0,n) Here n+1 means that the first one is init state. Cause here we need recursive calculate the cn based on previous cn, so add a init state for simplicity.
    cn(1) = 1;
    for k = 1:n
        for i = 1:N
            next = getNextd(RW(1:k,:,i),d);% Modified getNext for d dimension
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
mudBound = 2*d - 1 - 1/(2*d) - 3/(2*d)^2 - 16/(2*d)^3;
% mudu = 2d-1;
% mudl = d;

% Asymptotic bound on mud for large d found in Graham (2014); For d = 3, 4.6759, close to estimation.
% Result:
Ad = [1.21250812147076,1.14919779658453,1.22860321177269,1.23457689857401,1.14827765869649,1.19488297070291,1.32760060690718,1.17547408352294,1.29140037045228,1.21584621004824,1.18466070344711,1.28000952795186,1.21228964079189,1.22705964422744,1.17095313575684,1.19627036817860,1.28681299007487,1.32284826953376,1.19534287065679,1.16068704870388]
mud = [4.68436307067489,4.68300961627268,4.68445653100413,4.68251760043898,4.68442480157643,4.68319182266672,4.68326828837133,4.68209146451177,4.69074518532928,4.68392551232568,4.68356201712379,4.68250337027341,4.68140077136314,4.68464584198320,4.68300512887196,4.68307985938800,4.69138591576883,4.68893957447846,4.68191492431498,4.68589682718870]
gammad = [1.16660980980800,1.18270274096090,1.15357518348505,1.14990959855685,1.19088019024424,1.16431225504716,1.11834774338248,1.17760522571050,1.11951536786425,1.15642561470567,1.16216529760539,1.14043089517892,1.16667884447942,1.14567302386051,1.16987270513829,1.15934151988693,1.12084230101288,1.10465566844078,1.16681501835484,1.18414321917910]
%% 
%% problem 9. d dimension. SISR as problem 5, coefficients as problem 6. Similarly, coefficients estimation results are in bottom annotation. It is super slow...
big = 20;
slice = 3:big;
A = zeros(big,1);
mu = zeros(big,1);
gamma = zeros(big,1);
muBound = zeros(big,1);
muu = zeros(big,1);
mul = zeros(big,1);
for d = 3:20
%     d = 3; % Dimension. Be careful with d = 4, which is not supported here, since it has another Cn(d) expression.
    T = 1; % Estimation times
    Ad = zeros(T,1);
    mud = zeros(T,1);
    gammad = zeros(T,1);
    for t = 1:1 % 20 times estimation
        N = 2000;% Number of simulated particles. Reduced to 2000 since out of memory.
        % part = zeros(N,2); % init all particles, started at (0,0) in Z2 dimension. Column 1 represents for the x axis, and column 2 for y axis.
        n = 200;
        w = zeros(n+1,N); % init w. 
        w(1,:) = ones(1,N);
        RW = zeros(n+1,d,N); % Use multidimensional array to store generated Random Walks. And init all particles, started at (0,0) in Z2 dimension. Column 1 represents for the x axis, and column 2 for y axis, pages for different particles.
        cn = zeros(n+1,1); % Init cn. It is the value we want. Alternative: usedouble.empty(0,n) Here n+1 means that the first one is init state. Cause here we need recursive calculate the cn based on previous cn, so add a init state for simplicity.
        cn(1) = 1;
        for k = 1:n
            for i = 1:N
                next = getNextd(RW(1:k,:,i),d);% Modified getNext for d dimension
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
        if d == 4
%             x = 1:n;
%             X = x';
%             y = log(cn(2:end)) - 1/4.*log(x);
%             md = fitlm(X, y, 'Intercept', true); %Reference: https://www.mathworks.com/help/stats/linear-regression-model-workflow.html
%             coef = md.Coefficients.Estimate;
%             Ad(t) = exp(coef(1));
%             mud(t) = exp(coef(2));
%             gammad(t) = 0;
            Ad(t) = 0;
            mud(t) = 0;
            gammad(t) = 0;
        else
            x = 1:n;
            X = [x',log(x)'];
            y = log(cn(2:end));
            md = fitlm(X, y, 'Intercept', true); %Reference: https://www.mathworks.com/help/stats/linear-regression-model-workflow.html
            coef = md.Coefficients.Estimate;
            Ad(t) = exp(coef(1));
            mud(t) = exp(coef(2));
            gammad(t) = coef(3) + 1;
        end
    end
    A(d) = Ad(1);
    mu(d) = mud(1);
    gamma(d) = gammad(1);
    muBound(d) = 2*d - 1 - 1/(2*d) - 3/(2*d)^2 - 16/(2*d)^3;
    muu(d) = 2*d-1;
    mul(d) = d;
end
index = [3,5:big];
plot(index,mu(index),'r',index,muBound(index),'g',index,muu(index),'b--',index,mul(index),'b--');
legend('\mu_d','Asymptotic bound','2d-1','d')
% ax = gca;
set(gcf,'unit','centimeters','position',[10 5 20 12])
xlabel("d")
ylabel("\mu_d")
title("\mu_d estimation")

plot(index,A(index),'k',index,ones(1,length(index)))
legend("A_d","1")
% ax.XTickLabel(index)
set(gcf,'unit','centimeters','position',[10 5 20 12])
xlabel("d")
ylabel("A_d")
ylim([0 1.5])
title("A_d estimation")

plot(index,gamma(index),'k')
set(gcf,'unit','centimeters','position',[10 5 20 12])
xlabel("d")
ylabel("\gamma_d")
ylim([0 1.5])
title("\gamma_d estimation")