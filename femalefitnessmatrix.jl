function femalefitnessmatrix(m,n,s,hD,hN,σ)
    # Outputs the matrix Wf, where Wf[i,j] is the fitness of a female with genotype i,j

    na = (n-1)^m + 1;
    #print("na = ")
    #print(na)
    W = zeros(na);
    Wf = zeros(na,na);

    indexarray = linindices(m,n);
    #print("indexarray = ")
    #print(indexarray)

    nW = zeros(Int, na-1);
    nR = zeros(Int, na-1); #or na??

    for k = 1:na-1
        #print("haplotype = ")
        haplotypevector = indexarray[k,:];
        #print(haplotypevector)

        nn = zeros(Int, n);#changed from n-1

        for kk = 1:m
            nn[haplotypevector[kk]] +=1;
        end
        #print("number of alleles = ")
        #print(nn)
        if nn[3] == 0
            nW[k] = nn[1];
            #print("nW[k] = ")
            #print(nW[k])
            nR[k] = nn[2];
            #print("nR[k] = ")
            #print(nR[k])
        else nR[k] = -1;

        end

        w = 1;

        if sum(nn[3:end]) == 0 #no deleterious alleles present
            for kk = 1:m
                if haplotypevector[kk] == 2
                    w = w * (1 - σ[kk]);
                end
            end

        else w = 1-s;
        end

        W[k] = w;
        #print("W[k] = ")
        #print(W[k])
    end

    W[na] = 1-s;

    for i = 1:na
        for j = 1:na
            #print("i = ")
            #print(i)
            #print("j = ")
            #print(j)
            deleterious = zeros(Int, 2);
            if i == na || nR[i] == -1
                #print("yes for i")
                deleterious[1] = 1;
            end
            if j == na || nR[j] == -1
                #print("yes for j")
                deleterious[2] += 1;
                #print("deleterious[2] = ")
                #print(deleterious[2])
            end
            #print("deleterious?")
            #print(deleterious)

# QUESTION: Do the effects of multiple N alleles across target sites compound?

            if deleterious[1] == 1 && deleterious[2] == 1
                Wf[i,j] = 1-s;

            elseif i == na
                #before, this was nR[j], without "m - ". It had hD rather than hN.
                #The idea was that a single R allele confers partial resistance.
                #Now, all have to be R.
                if m - nR[j] == 0 #First heterozygous R...R/D...D
                    Wf[i,j] = (1 - hN*s) * W[j];
                elseif m - nR[j] > 0 #If any W alleles present, full cost
                                     #imposed
                    Wf[i,j] = (1 - hD*s) * W[j];
                end

            elseif j == na

                if m - nR[i] == 0
                    Wf[i,j] = (1 - hN*s) * W[i];
                elseif m - nR[i] > 0
                    Wf[i,j] = (1 - hD*s) * W[i]
                end

            elseif nR[i] == -1
                Wf[i,j] = (1 - hN*s) * W[j];
            elseif nR[j] == -1
                Wf[i,j] = (1 - hN*s) * W[i];
            else Wf[i,j] = W[i]*W[j];
            end

        end

    end

    return Wf

end
