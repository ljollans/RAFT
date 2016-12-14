function X = shuffle(A)  
tmp=randperm(length(A));
X=A(tmp);
end