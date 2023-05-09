function drivematrix(m,n,ϵ,ν,β)
    # Prints out the matrix κ. I+κ[i,j] is the probability of transition j->i under drive (ie excluding de novo mutations).

    na = (n-1)^m + 1;

    #kappasingle[i,j] is the probability of transition j->i during gametogenesis
    kappasingle = [1-ϵ 0.0 0.0 0.0;ϵ*ν*β 1.0 0.0 0.0;ϵ*ν*(1-β) 0.0 1.0 0.0;ϵ*(1-ν) 0.0 0.0 1.0];

    #κ[i,j] is the probability of transition j->i during gametogenesis (multilpex)
    κ = zeros(na,na)

    indexarray = linindices(m,n);

    for j=1:na-1

        #Introduce variable to keep track of total column probability, for use
        #later
        probtracker = 0.0;
        for i=1:na-1
            for kk=1:m
                #Calculate alleles at site m=kk for haplotypes i,j
                #print("ii = ")
                ii = indexarray[i,kk];
                #print(ii)
                #print("jj = ")
                jj = indexarray[j,kk];
                #print(jj)

                #First, the m=1 target site allele must transition j->i
                if kk==1
                    #print("kk = 1")
                    κ[i,j] = kappasingle[ii,jj];
                    #print("κ")
                    #print(i)
                    #print(j)
                    #print("=")
                    #print(κ[i,j])
                #Next, each subsequent allele must transition as well
                else
                    #print("kk =")
                    #print(kk)
                    κ[i,j] = κ[i,j]*kappasingle[ii,jj];
                    #print("κ")
                    #print(i)
                    #print(j)
                    #print("=")
                    #print(κ[i,j])
                end

            end

            probtracker = probtracker + κ[i,j];

        end

        #Finally, P(j->drive)=1-sum_{i<drive}P(j->i)
        κ[na,j] = 1 - probtracker;
    end

    #P(drive -> drive) = 1
    κ[na,na] = 1.0;

    #Subtract identity matrix
    for i=1:na
        κ[i,i] = κ[i,i] - 1;
    end

    return κ

end

#Checked!
