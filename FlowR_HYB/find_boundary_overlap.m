function boundary_list=find_boundary_overlap (A,partition_list,n)
list=zeros(n,1);
[row,col]=find(A);
indicator=partition_list(row)-partition_list(col);
[b_ind_row,~]=find(indicator);
boundary1=row(b_ind_row);
boundary2=col(b_ind_row);
list([boundary1;boundary2])=1;
boundary_list=find(list);
b_percentage=length(boundary_list)/n;
disp('boundary_nodes')
disp(length(boundary_list))
disp('percentage')
disp(b_percentage)
end
