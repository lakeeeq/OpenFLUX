function MDV = cauchy(a,b)

length_a = numel(a);
length_b = numel(b);

newDim = length_a + length_b - 1;

A = zeros(newDim);
B = [b';zeros(newDim-length_b,1)];
length_a_sub = length_a-1;
for i = 1:length_b
    A(i:length_a_sub+i,i) = a;
end
cc = 1;
for j = i+1:newDim
    A(j:length_a_sub+j-cc,j) = a(1:end-cc);
    cc = cc + 1;
end
MDV = A*B;
MDV = MDV';
end