function [E_extra,L_inv_cut,U_inv_cut,L_inv_com,U_inv_com,H4,H3,total_nodes,total_in,total_boundary,total_no_dup_b,total_no_dup] =Preprocess_overlap (partition_list,nparts,n,A,d,c)
if(~d)
    A=A+A';
else
    A=A';
end
vec = sum(A, 1);
vec = bsxfun(@max, vec, 1);
vec = 1 ./ vec;
D = spdiags(vec(:),0,n,n);
E=c*A*D;


%the boundary detection treats the matrix as undirected 
boundary_list=find_boundary_overlap (A,partition_list,n); %return the original indices of boundary nodes
%duplicating...
[community_package,dupli]=dup(boundary_list,partition_list,A,nparts,n,d);%return the boundary nodes for each cluster
L_inv=cell(nparts,1);
U_inv=cell(nparts,1);
L_inv_c=cell(nparts,1);
U_inv_c=cell(nparts,1);
wi=cell(nparts,1);
total_boundary=cell(1,nparts);
total_b_inner=cell(1,nparts);
total_nodes=cell(1,nparts);
total_in=cell(1,nparts);
total_no_dup=cell(1,nparts);
H3=cell(nparts,1);
H4=cell(nparts,1);
num_b_inner=zeros(nparts,1);
num_b=zeros(nparts,1);
num_total=zeros(nparts,1);
cluster_b_start=zeros(nparts,1);
dup_list=zeros(size(dupli));
member_bound_dup=cell(nparts,1);
for i=1:nparts
    member_in_i=find(community_package(i,:)==1);
    member_inner_bound_i=find(community_package(i,:)==2);
    member_bound_ori=find(community_package(i,:)==3);
    dup_i=partition_list(dupli)==i;
    if (sum(dup_i))
        [~,dup_loc]=ismember(dupli(dup_i),member_bound_ori);
        dup_list(dup_i)=dup_loc+length(member_inner_bound_i);
    end
    member_bound_dup{i}=find(community_package(i,:)==4);
    total_boundary{i}=[member_bound_ori,member_bound_dup{i}];
    total_b_inner{i}=[member_inner_bound_i,member_bound_ori];
    total_no_dup{i}=[member_in_i,total_b_inner{i}];
    total_nodes{i}=[total_no_dup{i},member_bound_dup{i}];
    I=speye(length(total_nodes{i}));
    wi{i}=I-E(total_nodes{i},total_nodes{i});
    [L_inv{i},U_inv{i}]=lu(wi{i}); 
    L_inv{i}=L_inv{i}^(-1); 
    U_inv{i}=U_inv{i}^(-1); 
    num_inner=length(member_in_i);
    num_b(i)=length(total_boundary{i});
    num_dup=length(member_bound_dup{i});
    num_b_inner(i)=length(total_b_inner{i});
    num_total(i)=length(total_nodes{i});
    cluster_b_start(i)=num_inner+length(member_inner_bound_i)+1;
    U_inv{i}=U_inv{i}(1:end-num_dup,:);
    if(num_b_inner(i))
        if (num_b(i)<size(L_inv{i},2))
            H3{i}=U_inv{i}(num_inner+1:end,:)*L_inv{i}(:,1:end-num_b(i));
            total_in{i}=total_nodes{i}(1:cluster_b_start(i)-1);
        else
            H3{i}=sparse(num_b_inner(i),1);
            total_in{i}=1;
        end
        H4{i}=U_inv{i}(num_inner+1:end,:)*L_inv{i}(:,end-num_b(i)+1:end);
    end
    if(num_b(i))
        L_inv_c{i}=L_inv{i}(:,cluster_b_start(i):end);
        U_inv_c{i}=U_inv{i};
        
    end
end
total_boundary=cell2mat(total_boundary);
total_b_inner=cell2mat(total_b_inner);
total_nodes=cell2mat(total_nodes);
total_in=cell2mat(total_in);
total_no_dup_b=cell2mat(total_no_dup(num_b~=0));
total_no_dup=cell2mat(total_no_dup);

actotal=cumsum(num_total);
biggest=actotal(end);
n_division=ceil(biggest/100000)+1;
n_i=zeros(n_division,1);
for i=2:n_division
    n_i(i)=find(actotal<=biggest/n_division*(i-1),1,'last');
end

L_inv_cut=cell(n_division,1);
U_inv_cut=cell(n_division,1);
U_inv_com=cell(n_division,1);
L_inv_com=cell(n_division,1);
for i=1:n_division-1
    L_inv_cut{i}=blkdiag(L_inv_c{n_i(i)+1:n_i(i+1)});
    U_inv_cut{i}=blkdiag(U_inv_c{n_i(i)+1:n_i(i+1)});
    L_inv_com{i}=blkdiag(L_inv{n_i(i)+1:n_i(i+1)});
    U_inv_com{i}=blkdiag(U_inv{n_i(i)+1:n_i(i+1)});
end
i=i+1;
L_inv_cut{i}=blkdiag(L_inv_c{n_i(i)+1:end});
U_inv_cut{i}=blkdiag(U_inv_c{n_i(i)+1:end});
L_inv_com{i}=blkdiag(L_inv{n_i(i)+1:end});
U_inv_com{i}=blkdiag(U_inv{n_i(i)+1:end});

b_inner_start=cumsum([1;num_b_inner]);
b_inner_end=b_inner_start(2:end)-1;
b_inner_start(end)=[];

b_start=cumsum([1;num_b]);
b_end=b_start(2:end)-1;
b_start(end)=[];

E_extra=E(total_boundary,total_b_inner);
for i=1:nparts
    E_extra(b_start(i):b_end(i),b_inner_start(i):b_inner_end(i))=0;
    if (~isempty(member_bound_dup{i}))
        group_i=partition_list(member_bound_dup{i});
        [~,dup_each_group]=ismember(member_bound_dup{i},dupli);
        pos=b_inner_start(group_i)+dup_list(dup_each_group)-1;
        E_extra(b_start(i):b_end(i),pos)=0;
    end
end
H3=blkdiag(H3{:});
H4=blkdiag(H4{:});
disp('Reduced boundary nodes')
disp(size(E_extra,1));
disp('percentage');
disp(size(E_extra,1)/n);
end

