% generate a symmetric Erdos-Renyi graph with given size N and link probability p
% output: adjacency matrix A and edge list E (.xls) 
function A = Generate_connected_Erdos_Renyi_graph(N,p)

A=zeros(N,N);

for i=1:(N-1)
    for j=(i+1):N
        if 0<=rand && rand<p
            A(i,j)=1;
            A(j,i)=1;
        end;
    end;
end;%first trial 

% while Connected_graph_check(A)==0
%     A=zeros(N,N);
%     for i=1:(N-1)
%         for j=(i+1):N
%             if 0<=rand && rand<p
%                A(i,j)=1;
%                A(j,i)=1;
%             end;
%         end;
%     end;
% end;% generate a connected graph A

% dim=0;
% for i=1:(N-1)
%     for j=(i+1):N
%         if A(i,j)==1
%             dim=dim+1;
%         end;
%     end;
% end;
% 
% E=zeros(dim,2);
% 
% % now express the graph A by edge list E
% k=1;
% for i=1:(N-1)
%     for j=(i+1):N
%         if A(i,j)==1
%             E(k,1)=i;
%             E(k,2)=j;
%             k=k+1;
%         end;
%     end;
% end;
% 
% xlswrite('graph_edge_list',E);


    

