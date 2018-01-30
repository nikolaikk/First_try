

a = rand(10000);
b = rand(10000);

tic

c = a * b;
c = c - mean(mean(c));

toc

imshow(c)