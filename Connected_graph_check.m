% Check if a given adjacency matrix corresponds to a connected graph
% Input: matrix A
% Output: 1 for connected; 0 for not connected

function re = Connected_graph_check(A)

N=length(A(1,:));
re=1;

% Path=zeros(N,N);
% 
% for i=1:N
%     for j=1:N
%         for k=1:N
%             Ak=A^k;
%             if Ak(i,j)>0
%                 Path(i,j)=1;
%                 break;
%             end;
%         end;
%     end;
% end;
% 
% for i=1:N
%     for j=1:N
%         if Path(i,j)==0
%             re=0;
%         end;
%     end;
% end;

sum=zeros(N,N);
for k=1:N
    sum=sum+A^k;
end;
for i=1:(N-1)
    for j=(i+1):N
        if sum(i,j)==0
            re=0;
            break;
        end;
        if sum(j,i)==0
            re=0;
            break;
        end;
    end;
end;



