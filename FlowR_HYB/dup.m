%execute after finding boundary, eliminate all nodes with cross edges of degree 1
%return community_package;
%each element of community_package is an array of length n;
%0 represents not belongs to this community
%1 represents inner nodes
%2 represents inner_boundary nodes
%3 represents non-duplicate boundary nodes
%4 represents duplicate boundary nodes
%delete degree one node, if only two nodes are connected,just duplicate one
function [community_package,dupli]=dup(boundary_list,partition_list,A,nparts,n,d)
if(d)
    A=A+A';
end
community_package=zeros(nparts,n);
dupli=zeros(n,1);
[I,reorder_b]=sort(partition_list(boundary_list));%reorder boundary nodes
indicator=diff([0;I]);
bound_start=find(indicator);
A_b_reorder=A(boundary_list(reorder_b),boundary_list(reorder_b));%get the reordered boundary node adj
l=length(boundary_list);
bound_end=[bound_start(2:end)-1;l];
j=1;
for i=1:length(bound_end)
    A_b_reorder(j:bound_end(i)-bound_start(i)+j,j:bound_end(i)-bound_start(i)+j)=0;
    j=bound_end(i)-bound_start(i)+j+1;
end
j=1;
for i=1:nparts
    community_package(i,partition_list==i)=1;
    b_list=partition_list(boundary_list)==i;
    if (sum(b_list))
        community_package(i,boundary_list(b_list))=3;
        A_b_reorder_i=A_b_reorder(bound_start(j):bound_end(j),:);%choose each cluster
        j=j+1;
        k=3;
        degree_list_i=sum(A_b_reorder_i,2)==1;
        if (sum(degree_list_i))
            A_b_reorder_one_i=A_b_reorder_i(degree_list_i,:);
            duplicate_list=sum(A_b_reorder_one_i,1)>1;
            if (sum(duplicate_list))
                duplicate_list_old_ordering=boundary_list(reorder_b(duplicate_list));
                dupli(duplicate_list_old_ordering)=1;
                community_package(i,duplicate_list_old_ordering)=4;
                [inner_bound_r,~]=find(A_b_reorder_one_i(:,duplicate_list));
                Ind_d=find(degree_list_i);
                inner_bound_r=boundary_list(reorder_b(Ind_d(inner_bound_r)+bound_start(j-1)-1));
                community_package(i,inner_bound_r)=2;
                
            end
        end
    end
end
dupli=find(dupli==1);
for i=1:nparts
    dup_ind=community_package(i,dupli)==3;
    if (sum(dup_ind))
        duped_nodes=dupli(dup_ind);
        nodes_i=find(community_package(i,:)==1);
        A_i=A(nodes_i,duped_nodes);
        [Ind,~]=find(A_i);
        community_package(i,nodes_i(Ind))=2;
    end
end
end




