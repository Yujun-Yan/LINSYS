function x = query_overlap (E_extra,L_inv_cut,U_inv_cut,L_inv_com,U_inv_com,L_k_inv,U_k_inv,T1,H3,total_nodes,total_in,total_boundary,total_no_dup_b,total_no_dup,s,n,c)
y=sparse(n,1);
x=zeros(n,1);
y(s)=1-c;
y_i=y(total_in);
y_b=y(total_boundary);
q=E_extra*(H3*y_i)+T1*y_b;
f=U_k_inv*(L_k_inv*q);
cuts=size(L_inv_cut,1);
e=cell(cuts,1);
ori=cell(cuts,1);
y_r=y(total_nodes);
j=1;
k=1;
for i=1:cuts
    l=size(L_inv_com{i},2);
    l2=size(L_inv_cut{i},2);
    ori{i}=U_inv_com{i}*(L_inv_com{i}*y_r(j:j+l-1));
    j=j+l;
    e{i}=U_inv_cut{i}*(L_inv_cut{i}*f(k:k+l2-1));
    k=k+l2;
end
x(total_no_dup)=cell2mat(ori);
x(total_no_dup_b)=cell2mat(e)+x(total_no_dup_b);
x=sparse(x);
end
