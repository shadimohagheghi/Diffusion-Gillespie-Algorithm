% This is the model in which social conversion happens before the self
% conversion

% two products

N=100;
p=0.4;

A = Generate_connected_Erdos_Renyi_graph(N,p);

while Connected_graph_check(A)==0
    A = Generate_connected_Erdos_Renyi_graph(N,p);
end;

A = [0,0,0,0,1;
     0,0,0,0,1;
     0,0,0,0,1;
     0,0,0,0,1;
     1,1,1,1,0];
% generate a N*N complete graph
for i=1:N
    for j=1:N
        if i~=j
            A(i,j) = 1;
        end;
    end;
end;

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
alpha(N) = 0.5;

Delta = zeros(2,2);
Delta(1,1) = rand; Delta(1,2) = 1-Delta(1,1);
Delta(2,1) = rand; Delta(2,2) = 1-Delta(2,1);

Delta = [0, 1;
         1, 0];

% define P
T=500;
P=zeros(N,T);

% generate the initial condition randomly
P(:,1) = rand(N,1);
%P(:,1) = [0;0;0;0.01;0];
P(:,1) = zeros(N,1);
P(N,1) = 1;

for t=1:(T-1)
    P(:,t+1) = Delta(2,2)*( eye(N)-diag(alpha) )*P(:,t) + ( Delta(1,2)*eye(N) + Delta(1,1)*diag(alpha) )*A_tilde*P(:,t) + ( Delta(2,1)-Delta(1,2) )*( eye(N) - diag(alpha) )*diag( P(:,t) )*A_tilde*P(:,t);
end;


