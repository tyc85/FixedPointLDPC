load('G_trim_and_H_indx.mat');

u = rand(1, 1978) < 0.5;

p = mod(u*G_Trim, 2); %generate parities

c = zeros(2209, 1);
% need to assign to the right order
for i = 1:1978
	c(indx_info(i)) = u(i);
end


for i = 1:231
	for j = 1:1978
		c(indx_par(i));
	end
end

mod(H*c, 2)