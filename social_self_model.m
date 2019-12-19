% social-self model
% compare the solution to the approximation model and the Markov-chain
% model

N = 10;
R = 4;
%R = 2;
p = 0.5;

T = 50;

% generate the graph
A = Generate_connected_Erdos_Renyi_graph(N,p);

while Connected_graph_check(A)==0
    A = Generate_connected_Erdos_Renyi_graph(N,p);
end;

% % generate a scale-free network
% seed =[0 1 0 0 1;1 0 0 1 0;0 0 0 1 0;0 1 1 0 0;1 0 0 0 0];
% A = SFNG(N, 3, seed);


% generate A manually
% A = [0,0,0,0,1;
%      0,0,0,0,1;
%      0,0,0,0,1;
%      0,0,0,0,1;
%      1,1,1,1,0];

A = [0,0,0,0,0,0,0,0,0,1;
     0,0,0,0,0,0,0,0,0,1;
     0,0,0,0,0,0,0,0,0,1;
     0,0,0,0,0,0,0,0,0,1;
     0,0,0,0,0,0,0,0,0,1;
     0,0,0,0,0,0,0,0,0,1;
     0,0,0,0,0,0,0,0,0,1;
     0,0,0,0,0,0,0,0,0,1;
     0,0,0,0,0,0,0,0,0,1;
     1,1,1,1,1,1,1,1,1,0];
    

% generate A_tilde from A
Row_sum = zeros(1,N);
for i=1:N
    Row_sum(i) = sum( A(i,:) );
end;

A_tilde = zeros(N,N);

for i=1:N
    A_tilde(i,:) = A(i,:)/Row_sum(i);
end;

% generate the model parameters
 alpha = rand(N,1);
% manually set alpha
% alpha = [0.1;0.2;0.3;0.4;0.5];
% alpha = 0.95*ones(N,1);

Delta = rand(R,R);
for r = 1:R
    Delta(r,:) = Delta(r,:)/sum(Delta(r,:));
    Delta(r,R) = 1;
    for s = 1:(R-1)
        Delta(r,R) = Delta(r,R)-Delta(r,s);
    end;
end;

Delta = [0.6,0.4,0,0;
         0.3,0.7,0,0;
         0,0,1,0;
         0,0.8,0,0.2];

% % manually set Delta
% Delta = [0.2,0.3,0.4,0.1;
%          0.6,0.1,0.1,0.2;
%          0.3,0.3,0.1,0.3;
%          0.2,0.1,0.6,0.1];



% Part 1: Approximation solution
% P(:,:,t): trajectory matrix

% define the trajectory
P = zeros(N,R,T);
% generate the initial condition
P(:,:,1) = rand(N,R);
for i = 1:N
    P(i,:,1) = P(i,:,1)/sum(P(i,:,1));
    P(i,R,1) = 1;
    for s = 1:(R-1)
        P(i,R,1) = P(i,R,1) - P(i,s,1);
    end;
end;

% start iteration
for t = 1:(T-1)
    P(:,:,t+1) = diag(alpha)*A_tilde*P(:,:,t) + ( eye(N) - diag(alpha) )*P(:,:,t)*Delta;
    for i = 1:N
        P(i,R,t+1) = 1;
        for s = 1:(R-1)
            P(i,R,t+1) = P(i,R,t+1) - P(i,s,t+1);
        end;
    end;
end;

Q = zeros(N,R,T);% count the # in K trials
K = 100000;
for k=1:K
    
    s = zeros(N,1); % vector D(t+1)
    s_pre = zeros(N,1); % vector D(t)
    
    for i = 1:N % 4 possible states
        % first generate the initial condition D(1)
        a = rand;
        if 0 <= a && a < P(i,1,1)
            s(i) = 1;
            Q(i,1,1) = Q(i,1,1) + 1;
            continue;
        end;
        
        if P(i,1,1) <= a && a < P(i,1,1)+P(i,2,1)
            s(i) = 2;
            Q(i,2,1) = Q(i,2,1) + 1;
            continue;
        end;
        
        if P(i,1,1)+P(i,2,1) <= a && a < 1-P(i,4,1)
            s(i) = 3;
            Q(i,3,1) = Q(i,3,1) + 1;
            continue;
        end;
        
        if 1-P(i,4,1) <= a && a < 1
            s(i) = 4;
            Q(i,4,1) = Q(i,4,1) + 1;
            continue;
        end;
    end;
    
    
    
    
        for t=2:T
            s_pre = s;
            for i = 1:N
                % first social conversion
                % randomly pick a neighbor j
                j = unidrnd(N);
                while A(i,j) == 0
                    j = unidrnd(N);
                end;
                % follow j's previous state with prob. alpha_i
                social = rand;
                if 0 <= social && social < alpha(i)
                    s(i) = s_pre(j);
                    Q(i,s(i),t) = Q(i,s(i),t) + 1;
                    continue;
                end;
                
                if alpha(i) <= social && social < 1 % if not follow (with prob. 1-alpha_i)
                    self = rand; % start self conv
                    if 0 <= self && self < Delta(s(i),1)
                        s(i) = 1;
                        Q(i,1,t) = Q(i,1,t) + 1;
                        continue;
                    end;
                    if Delta(s(i),1) <= self && self < Delta(s(i),1) + Delta(s(i),2)
                        s(i) = 2;
                        Q(i,2,t) = Q(i,2,t) + 1;
                        continue;
                    end;
                    if Delta(s(i),1) + Delta(s(i),2) <= self && self < 1 - Delta(s(i),4)
                        s(i) = 3;
                        Q(i,3,t) = Q(i,3,t) + 1;
                        continue;
                    end;
                    if 1 - Delta(s(i),4) <= self && self < 1
                        s(i) = 4;
                        Q(i,4,t) = Q(i,4,t) + 1;
                        continue;
                    end;
                end;
            end;
        end;


% % This part is specially for the star network
%     for t=2:T
%         s_pre = s;
%         for i = 1:(N-1) % peripheral nodes
%             % social conversion, can only pick N
%             social = rand;
%             if 0 <= social && social < alpha(i)
%                s(i) = s_pre(N);
%                Q(i,s(i),t) = Q(i,s(i),t) + 1;
%                continue;
%             end;
%             if alpha(i) <= social && social < 1 % if not follow (with prob. 1-alpha_i)
%                self = rand; % start self conv
%                if 0 <= self && self < Delta(s(i),1)
%                   s(i) = 1;
%                   Q(i,1,t) = Q(i,1,t) + 1;
%                   continue;
%                end;
%                if Delta(s(i),1) <= self && self < Delta(s(i),1) + Delta(s(i),2)
%                   s(i) = 2;
%                   Q(i,2,t) = Q(i,2,t) + 1;
%                   continue;
%                end;
%                if Delta(s(i),1) + Delta(s(i),2) <= self && self < 1 - Delta(s(i),4)
%                   s(i) = 3;
%                   Q(i,3,t) = Q(i,3,t) + 1;
%                   continue;
%                end;
%                if 1 - Delta(s(i),4) <= self && self < 1
%                   s(i) = 4;
%                   Q(i,4,t) = Q(i,4,t) + 1;
%                continue;
%                end;                 
%             end;            
%         end;
%         % now for node n
%         j = unidrnd(N);
%         while A(i,j) == 0
%               j = unidrnd(N);
%         end;
%         social = rand;
%         if 0 <= social && social < alpha(N)
%            s(N) = s_pre(j);
%            Q(N,s(N),t) = Q(N,s(N),t) + 1;
%            continue;
%         end;
%         if alpha(N) <= social && social < 1 % if not follow (with prob. 1-alpha_i)
%                self = rand; % start self conv
%                if 0 <= self && self < Delta(s(N),1)
%                   s(N) = 1;
%                   Q(N,1,t) = Q(N,1,t) + 1;
%                   continue;
%                end;
%                if Delta(s(N),1) <= self && self < Delta(s(N),1) + Delta(s(N),2)
%                   s(N) = 2;
%                   Q(N,2,t) = Q(N,2,t) + 1;
%                   continue;
%                end;
%                if Delta(s(N),1) + Delta(s(N),2) <= self && self < 1 - Delta(s(N),4)
%                   s(N) = 3;
%                   Q(N,3,t) = Q(N,3,t) + 1;
%                   continue;
%                end;
%                if 1 - Delta(s(N),4) <= self && self < 1
%                   s(N) = 4;
%                   Q(N,4,t) = Q(N,4,t) + 1;
%                   continue;
%                end;                 
%             
%         end;     
%         
%     end;     
        
    k
end;
Q = Q/K;   

% alpha
% 
% figure;
% for i=1:N
% for r=1:R
% hold on;
% pir = zeros(T,1);
% for t=1:T
% pir(t) = Q(i,r,t);
% end;
% plot([1:T],pir);
% end;
% end;

p12 = zeros(T,1);
q12 = zeros(T,1);

for t = 1:T
    p12(t) = mean(P(:,2,t));
    q12(t) = mean(Q(:,2,t));
end;

figure;
plot([1:T],p12,'b');
hold on;
plot([1:T],q12,'r');
                   
                    



