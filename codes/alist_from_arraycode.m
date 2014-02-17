p = 47;
r = 5;
I = eye(p);
S = circshift(I,1);

H = kron(ones(1, p), I);
row = eye(p);
for i = 2:r
	for j = 2:p
		row = [row S^(i*(j -1))];
	end
	H = [H;row];
	row = eye(p);
end
% a_binary=gf(H); % Turn 'a' into a special binary variable.
% rank(a_binary)

FID = fopen('H_array_p47_r5.txt', 'w');
fprintf(FID, '%d %d\n', p^2, r*p);
fprintf(FID, '%d %d\n', r, p);
for i = 1:p^2
	fprintf(FID, '%d ', r);
end
fprintf(FID, '\n');
for i = 1:p*r
	fprintf(FID, '%d ', p);
end

fprintf(FID, '\n');
for i=1:p^2
	ind = find(H(:, i));
	for j = 1:r
		fprintf(FID, '%d ', ind(j)- 1);
	end
	fprintf(FID, '\n');
end

for i=1:p*r
	ind = find(H(i, :));
	for j = 1:p
		fprintf(FID, '%d ', ind(j) - 1);
	end
	fprintf(FID, '\n');
end
	
fclose(FID);
	
