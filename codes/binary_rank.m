a=[1 1 0 1 1 0 0;...
   0 1 1 0 1 0 0;...
   0 0 1 1 1 1 1;...
   1 0 0 0 1 1 1];
a_binary=gf(a); % Turn 'a' into a special binary variable.
rank(a_binary) % Find the rank using binary operations.