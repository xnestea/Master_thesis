%module add matlab/R2018b
%matlab -nodisplay -nodesktop
%change threshold and name of input / output files
%cd /home/xnestea/Master/Data/REO_to_K/Identify_stable_gene_pair/Binomial_test
load ("merged_data_prep.mat")

cb1_res = sig_stable_gene_pair(rownames,T1.data,0.0000000001);
T_rest = [T2.data T3.data];
ca_23_res = sig_stable_gene_pair(rownames,T_rest, 0.0000000001);
c1_res = c1_res{1,1};
c_23_res = c_23_res{1,1};
c_23_res(:,[1 2]) = c_23_res(:,[2 1]);
C1 = intersect(c1_res, c_23_res, "rows");
clear c1_res c_23_res T_rest

c2_res = sig_stable_gene_pair(rownames,T2.data,0.0000000001);
T_rest = [T1.data T3.data];
c_13_res = sig_stable_gene_pair(rownames, T_rest,0.0000000001);
c2_res = c2_res{1,1};
c_13_res = c_13_res{1,1};
c_13_res(:,[1 2]) = c_13_res(:,[2 1]);
C2 = intersect(c2_res, c_13_res, "rows");
clear c2_res c_13_res T_rest

c3_res = sig_stable_gene_pair(rownames,T3.data,0.0000000001);
T_rest = [T1.data T2.data];
c_12_res = sig_stable_gene_pair(rownames, T_rest,0.0000000001);
c3_res = c3_res{1,1};
c_12_res = c_12_res{1,1};
c_12_res(:,[1 2]) = c_12_res(:,[2 1]);
C3 = intersect(c3_res, c_12_res, "rows");
clear T_rest T1 T2 T3 colnames1 colnames2 colnames3 rownames c3_res c_12_res

save("Merged_7GP.mat")