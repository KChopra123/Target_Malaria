function linindices(m,n)
    # Returns the ((n-1)^m+1)xm array [1 1 1 ... 1;2 1 ... 1;3 1 ... 1; ... ; n-1 1 ... 1;2 1 ... 1;2 2 1 ... 1; ... ;n-1 n-1 ... n-1; n n n]
    # Try, say, linindices(3,4) 

    na = (n-1)^m+1;

    linindicesarray = zeros(Int,na,m);

    for i = 1:na-1
        for j = 1:m
            linindicesarrayindex = ceil(Int, i/((n-1)^(j-1)));
            linindicesarray[i,j] = mod(linindicesarrayindex, (n-1));

            if linindicesarray[i,j] == 0
                linindicesarray[i,j] = n-1;
            end

        end
    end

    for j = 1:m
        linindicesarray[na,j] = n;
    end

    return linindicesarray
end
