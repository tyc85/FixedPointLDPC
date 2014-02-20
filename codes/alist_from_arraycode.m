p = 47;
r = 5;
rate = 1 - (r*p - r + 1)/p^2
I = eye(p);
% shift backward
% S = circshift(I,1);
% shift forward
S = circshift(I,p-1);
H = kron(ones(1, p), I);
row = eye(p);
proto = zeros(5, 47);
for i = 2:r
	for j = 2:p
		proto(i, j) = (i-1)*(j -1);
	end
end
proto;
for i = 2:r
	for j = 2:p
		row = [row S^((i-1)*(j -1))];
	end
	H = [H;row];
	row = eye(p);
end
% a_binary=gf(H); % Turn 'a' into a special binary variable.
% rank(a_binary)

FID = fopen('H_array_p47_r5_forward.txt', 'w');
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

% output columns
for i=1:p^2
	ind = find(H(:, i));
	for j = 1:r
		fprintf(FID, '%d ', ind(j)- 1);
	end
	fprintf(FID, '\n');
end


% output rows
for i=1:p*r
	ind = find(H(i, :));
	for j = 1:p
		fprintf(FID, '%d ', ind(j) - 1);
	end
	fprintf(FID, '\n');
end
	
fclose(FID);
	
