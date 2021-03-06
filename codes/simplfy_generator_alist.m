
 p = 47;
 r = 5;
% I = eye(p);
% S = circshift(I,1);

% H = kron(ones(1, p), I);
% row = eye(p);
% for i = 2:r
% 	for j = 2:p
% 		row = [row S^(j -1)];
% 	end
% 	H = [H;row];
% 	row = eye(p);
% end
FID2 = fopen('column_flag_arrary_p47_r5.txt', 'r');
column_flag = fscanf(FID2, '%d');
% note the index with 0 gives the index in H for infobit
% it seems to confusing but it's correct=> index_par is the identity
% part after Gaussian elimination
indx_info = find(1- column_flag); % find zero
indx_par = find(column_flag); % find 1
fclose(FID2);
FID = fopen('H_array_p47_r5g.txt', 'r');

all = fscanf(FID, '%d');
n = all(1);
m  = all(2);
% n = psqr;
% m = pr;
vdeg_max = all(3);
cdeg_max = all(4);
counter = 5;
vdeg = zeros(1, n);
cdeg = zeros(1, m);

% p = sqrt(psqr);
% r = rp/p;
for i = 1:p^2
	vdeg(i) = all(counter);
	counter = counter + 1;
end

for i = 1:p*r
	cdeg(i) = all(counter);
	counter = counter + 1;
end

vlist = zeros(n, max(vdeg));
clist = zeros(m, max(cdeg));


for i=1:p^2
	for j = 1:vdeg(i)
		vlist(i, j) = all(counter);
		counter = counter + 1;
	end
end

for i=1:m
	for j = 1:cdeg(i)
		clist(i, j) = all(counter);
		counter = counter + 1;
	end
end
% vlist is the columns in the generator matrix that corresponds to parities
% the element in vlist is the row index where the entry is non-zero
% for i = 1:numel(indx_par)
% 	indx = indx_par(i);
% 	for j = 1:vdeg(indx)
% 		G_Trim(i, vlist(indx, j)) = 1;
% 	end
% end
% need to interleave the parity into right order 
fclose(FID);
%%
FID_out = fopen('G_array_p47_r5.txt', 'w');

fprintf(FID_out, '%d %d', n-m, m);%  # of information bits and # of parities
fprintf(FID_out, '\n');
% output the index for info bits, the index is reduced by 1 to match c++
% indexing
for i = 1: numel(indx_info)
	fprintf(FID_out, '%d ', indx_info(i)-1);
end
fprintf(FID_out, '\n');


% output the index for parity bits, the index is reduced by 1 to match c++
% indexing
for i = 1: numel(indx_par)
	fprintf(FID_out, '%d ', indx_par(i)-1);
end
fprintf(FID_out, '\n');

for i=1:m
	for j = 1:cdeg(i)
		fprintf(FID_out, '%d ', clist(i, j));
	end
	fprintf(FID_out, '\n');
end
	

fclose(FID_out);
	
