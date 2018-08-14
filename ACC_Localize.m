function L = ACC_Localize(n,r)

L = zeros(n,n);


for i=1:n
    for j=i+1:n
        dist = min(abs(i-j),mod(-abs(i-j),n));
        L(i,j) = exp(-dist*dist/2/r/r);
    end
end

L = eye(n) + L + L';

end

