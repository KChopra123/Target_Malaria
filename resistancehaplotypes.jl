function resistancehaplotypes(m,n,resistancealleles)
# if resistancealleles is a vector of the resistant allele positions,
# resistancehaplotypes calculates the resistant haplotype positions
    # eg. if R,N are resistant and m=2, RR,RN,NR,NN are resistant haplotypes

    na = (n-1)^m + 1;
    resistancehaplotypes = zeros(Int,na);

    indexarray = linindices(m,n);


    for iiii = 1:na
        checker = 1;
        for jjjj = 1:m
            if (indexarray[iiii,jjjj] in resistancealleles) == false
                checker = 0;
            end
        end

        if checker == 1
            resistancehaplotypes[iiii] = 1;
        end
    end

    return resistancehaplotypes

end
