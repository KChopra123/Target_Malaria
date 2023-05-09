function allele_comb(n)
    length = n^2;
    v = zeros(Int,length,2);

    position = 1;
    for i = 1:n
        for j = 1:n
            v[position,:] = [i,j];
            position += 1;
        end
    end

    return v
end

function genotype_position(k1,k2,n)
    return (k1-1)*n + k2
end

function mating(females,males,pmat)
    # Takes vectors of male & female genotypes, and produces vector of mated females
    # Given one male and one female, pmat is the probability that they find each other and mate (pmat does not account for fecundity)

    # Check number of male genotypes = number of female genotypes
    if size(males) != size(females)
        println("Error - males and umated female vectors should be of the same length")
        return
    end

    # Let n denote the number of genotypes
    n = size(males)[1]
    #println(n)

    tot_males = sum(males);
    #println("tot_males = $(tot_males)")
    #tot_females = sum(females);

    # One entry in the vector of mated females for each unordered pair of genotypes, and one for females which don't find a mate.
    # We count ordered pairs of genotypes in the following way, as given by the function allele_comb above.
    # Consider n genotypes called 1,2,3,...,n. Order pairs as (1,1),...,(1,n),(2,1),...,(2,n),...,(n,1),...,(n,n).
    # Including females which do not find a mate, there are (1/2)n(n+1)+1 entries. The final entry is non-mated females.
    mated_females = zeros(n^2 + 1);
    #println("size of mated females vector is $(size(mated_females))")

    # If there are no males, there will be no mated females
    if tot_males == 0
        return mated_females
    end

    # The vector prob_vec lists the probabilities of mating with a male of each genotype.
    # The entry prob_vec[i], for i≤n, denotes the probability of mating with a male of genotype i.
    # The entry prob_vec[n+1] denotes the probability of not finding a mate, given by (1 - pmat)^tot_males.
    prob_vec = zeros(n+1);
    for k = 1:n
        prob_vec[k] = (1 - (1 - pmat)^tot_males) * males[k]/tot_males;
    end
    prob_vec[end] = (1 - pmat)^tot_males;
    #println("Probability vector is $(prob_vec), with sum $(sum(prob_vec))")

    # Now find the (random) number of females with genotype 1≤i≤n which mate with males of genotype 1≤j≤n, or not at all.
    # This will gve a mated female with genotype pair (i,j).
    # The first entry in the pair corresponds to the female's genotype, the second to the male's.
    # Keeping track of which genotype originates from which sex is important, as only female reproductive fitness differs between genotypes.
    # This affects the number of eggs produced, not the probability of an individual surviving.
    # Unmated females in the entry mated_females[n+1] will obviously not reproduce.
    for i = 1:n
        v = rand(Multinomial(females[i],prob_vec));
        for j = 1:n
            mated_females[genotype_position(i,j,n)] += v[j];
        end
        mated_females[end] += v[end];
    end

    return mated_females

end

function genotype_position_unordered(k1,k2)
    if k1 < k2
        (k1,k2) = (k2,k1);
    end
    return trunc(Int, 0.5*k1*(k1-1) + k2)
end

function allele_comb_unordered(n)
    length = trunc(Int, 0.5*n*(n+1));
    v = zeros(Int,length,2);

    position = 1;
    for i = 1:n
        for j = 1:i
            v[position,:] = [i,j];
            position += 1;
        end
    end

    return v
end

function gamete_freq(genotype,ϵ,ν,μ,β,ξ)
    # The simulation runs for r=4 alleles (one gRNA). Individuals have one of n=(1/2)r(r+1) genotypes.
    # The program outputs the probability vector v.
    # Given a genotype, v[i] denotes the probability that a gamete produced has allele i.

    r = 4;
    n = trunc(Int, 0.5*r*(r+1));
    v = zeros(r);

    # The n-vector w_j contains the probabilities that the j^th allele of the genotype will produce a gamete of type i.
    # With Mendelian inheritance, and without mutation, w_j would trivially be 1 for the entry of the j^th allele of the genotype.
    # Mutation and drive make this non-tivial, eg. a wild allele can be converted in the prescence of drive, and hence produce non-wild gametes.
    w1 = zeros(r);
    w2 = zeros(r);

    n = 10; # Number of genotypes
    # Pairs of alleles within a genotype are unordered, ie. genotype (i,j) is equivalent to (j,i).
    # So we may consider (i,j) and (j,i) to be the same.
    # We order the genotypes as follows (as in the programs genotype_position_unordered, allele_comb_unordered)
    # (1,1), (2,1), (2,2), (3,1), (3,2), (3,3), (4,1),...
    # So eg. genotype 4 has alleles of type 1 and 3.
    # The vector alleles contains the pair of alleles in the genotype.
    # By construction, alleles[1] ≥ alleles[2].
    alleles = allele_comb_unordered(r)[genotype,:];
    if maximum(alleles) > 4
        print("Error - alleles must have entries from 1 ro 4")
        return -0.5
    end

    # First consider de novo mutations. Wild type (type 1) alleles spontaneously mutate:
    #  - to R (type 2) with probability ξ*μ
    #  - to N (type 3) with probability (1-ξ)*μ
    if alleles[1] < 4   # and hence alleles[2]≤alleles[1]<4, so no drive
        if alleles[2] > 1   # so no wild type
            v[alleles[1]] += 0.5;
            v[alleles[2]] += 0.5;
            return v
        elseif alleles[1] > 1   # So one wild type and one R or N
            w1[alleles[1]] = 1;
            w2 = [1-μ;ξ*μ;(1-ξ)*μ;0];
            v = 0.5*w1 + 0.5*w2;
            return v
        else    # Two wild type alleles
            v = [1-μ;ξ*μ;(1-ξ)*μ;0];
            return v
        end
    elseif alleles[2] == 1  # So one wild type, and one drive allele
    # Now have to take into account drive, and therefore NHEJ
    # See PNAS paper for more info on this
        w1[4] = 1;
        w2 = [(1-μ)*(1-ϵ);ξ*μ+(1-μ)*β*ϵ*ν;(1-ξ)*μ+(1-μ)*(1-β)*ϵ*ν;(1-μ)*ϵ*(1-ν)];
        v = 0.5*w1 + 0.5*w2;
        return v
    else    # So one wild type and one R or N. Inheritance is Mendelian.
        v[alleles[1]] += 0.5;
        v[alleles[2]] += 0.5;
        return v
    end

    print("Error with gamete_freq function - options exhausted.")

end

function gamete_freq_array(ϵ,ν,μ,β,ξ)
    # Returns a 10x4 array. The (i,j)^th entry is the probability that a gamete produced by an individual of genotype has allele j.
    r = 4;
    n = 10;
    gamete_array = zeros(n,r);
    for i = 1:n
        gamete_array[i,:] = transpose(gamete_freq(i,ϵ,ν,μ,β,ξ));
    end

    return gamete_array
end

function female_fitness(s,hD,hN,σ)
    # Outputs fitness_vector. fitness_vector[i] is the reproductive fitness of a female with genotype i.
    # This program really only works for one target site
    r = 4;
    n = 10;

    # Recall from above the ordering of genotypes: (1,1),(2,1),(2,2),(3,1),...
    allele_list = allele_comb_unordered(r);
    n = size(allele_list)[1];
    fitness_vector = zeros(n);

    for i = 1:n
        alleles = allele_list[i,:];
        if alleles[2] > 2    # Therefore both alleles are N or D
            fitness_vector[i] = 1-s;
        else
            heterozygote_fitness = [1;1-σ;1-hN;1-hD];
            fitness_vector[i] = heterozygote_fitness[alleles[1]] * heterozygote_fitness[alleles[2]];
        end
    end

    return fitness_vector
end

function offspring_production(mated_genotype_pair,no_females,totfem,Rm,K,fitness_vector,gamete_array)
    # Generates the genotypes of the offspring in generation t+1 from a pair of genotypes in a mated female - ie. her genotype and genotype of sperm donor.
    # Outputs two vectors - female_offspring and male_offspring
    # female_offspring[i] is the number of females of genotype i produced who reach adulthood. These may then mate. Similar for male_offspring.
    # individuals denotes the number of individuals of that genotype in generation t.
    # no_females = number of females in generation t, including all genotypes
    # totpop = total population in generation t.
    # Rm = max. reproduction rate. Corresponds to average number of offspring without competition
    # K = carrying capacity in Beverton-Holt model
    # fitness_vector[i] is the reproductive fitness of a female of genotype i.
    # If a female of fitness 1 produces M offspring in expectation, a female of fitness produces M*f offspring in expectation.
    # gamete_array[i,j] is the probability that a gamete produced by an individual of genotype i has allele/haplotype j.
    n = 10; # n denotes the number of genotypes

    # First ascertain genotypes of mother and father
    (female_genotype,male_genotype) = allele_comb(n)[mated_genotype_pair,:];    # Streamline this later if necessary

    # First calculate the *female* fitness
    fitness = fitness_vector[female_genotype];

    # Calculate the number of eggs produced - use Beverton-Holt model
    # If the female fitness is zero, no male or female offspring of any genotype are produced
    if fitness == 0
        return(zeros(Int,n),zeros(Int,n))
    else
        exptd_pop = 2*fitness*Rm*no_females/(1+2*totfem*(Rm-1)/K);
        #exptd_pop = fitness*Rm*no_females/(1+totfem*(Rm-1)/K);
        quantity_female_offspring = rand(Poisson(0.5*exptd_pop));
        quantity_male_offspring = rand(Poisson(0.5*exptd_pop));
    end

    # Now assign genotypes to the eggs
    # Let n denote the number of genotypes
    n = size(fitness_vector)[1];
    genotype_ordering = allele_comb_unordered(4);
    #female_alleles = genotype_ordering[female_genotype,:];
    #male_alleles = genotype_ordering[male_genotype,:];

    # Let offspring_prob[i] denote the probability that an individual offspring has genotype i.
    offspring_prob = zeros(n);
    for i = 1:n
        alleles = genotype_ordering[i,:];   # Alleles (of offspring) in genotype i
        if alleles[1] == alleles[2] # Eg. for offspring to be WW, but parental gametes must have W
            offspring_prob[i] = gamete_array[female_genotype,alleles[1]] * gamete_array[male_genotype,alleles[2]];
        else    # Eg. for offspring to be WR, either (W,R) or (R,W) for (mother's gamete,father's gamete)
            offspring_prob[i] = gamete_array[female_genotype,alleles[1]] * gamete_array[male_genotype,alleles[2]] + gamete_array[female_genotype,alleles[2]] * gamete_array[male_genotype,alleles[1]];
        end
    end

    female_offspring = rand(Multinomial(quantity_female_offspring,offspring_prob));
    male_offspring = rand(Multinomial(quantity_male_offspring,offspring_prob));

    return(female_offspring,male_offspring)

end

function TSR_mated_mixed(s,hD,hN,σ,Rm,K,xD,yD,ϵ,ν,μ,β,ξ,pmat,preexisting,T::Int,test_threshold)
    # Analogue of TargetSiteResistance program for the modified mated females regime. Used to verify algorithm against original.
    # Outputs (x,y,popsize,xnum,ynum,resistance,τ). x,y[i,t] = proportion of *allele* i at time t. xnum,ynum[i,t] = number of allele i. Other outputs as TargetSiteResistance.
    # For m=1 target site only, for now.
    # Don't have an estimate for preexisting variation, so will just have to run burn-in for a long time to obtain standing variation.

    #Set the number of alleles, r, and the number of genotypes, n = (1/2)r(r+1)
    r = 4;
    n = 10;
    # genotype_list[i,:] contains the alleles within genotype i.
    genotype_list = allele_comb_unordered(r);

    # First calculate gamete frequencies for the genotypes (see gamete_freq_array function for more details) and females fitnesses (see female_fitness funtion for details)
    gamete_array = gamete_freq_array(ϵ,ν,μ,β,ξ);
    fitness_vector = female_fitness(s,hD,hN,σ);

    # At first all individuals have genotype 1 (homozygous W)
    females = zeros(Int,n);
    males = zeros(Int,n);
    females[1] = rand(Poisson(0.5*K));
    males[1] = rand(Poisson(0.5*K));

    # Now burn-in if wanted
    if preexisting == 1
        for t = 1:ceil(Int,5/σ)   # Longer burn-in time than initial program as we don't start from an estimate of the standing variation
            mated_females = mating(females,males,pmat);
            # Throw away females which did not find a mate.
            mated_females = mated_females[1:end-1];

            # Migration does not occur in this algorithm
            # So skip to reproduction
            # Now for females and males in next generation
            totfem = sum(females);
            females = zeros(Int,n);
            males = zeros(Int,n);
            for kk = 1:n^2
                (female_genotype_offspring,male_genotype_offspring) = offspring_production(kk,mated_females[kk],totfem,Rm,K,fitness_vector,gamete_array);
                females += female_genotype_offspring;
                males += male_genotype_offspring;
            end
        end
    end

    # Main phase of simulation
    # Will output (x,y,popsize,xnum,ynum,resistance,τ)
    xnum = zeros(r,T);
    ynum = zeros(r,T);
    x = zeros(r,T);
    y = zeros(r,T);
    popsize = zeros(3,T);
    # Initialise resistance, τ
    resistance = 0;
    τ = NaN;

    popsize[:,1] = [sum(females);sum(males);sum(females+males)];
    # Each allele contributes 1 towards xnum. So an RR individual contributes 2 to xnum[2,t].
    # This means that sum(xnum[:,t]) = 2*sum(females). We will halve xnum, ynum at the end to correct for this.
    for i = 1:n
        xnum[genotype_list[i,1],1] += females[i];
        xnum[genotype_list[i,2],1] += females[i];
        ynum[genotype_list[i,1],1] = males[i];
        ynum[genotype_list[i,2],1] = males[i];
    end
    x[:,1] = xnum[:,1];
    y[:,1] = ynum[:,1];

    if maximum(x[:,1]) > 0
        x[:,1] = x[:,1]/sum(x[:,1]);
    end
    if maximum(y[:,1]) > 0
        y[:,1] = y[:,1]/sum(y[:,1]);
    end

    # Now add drive. For now, all drive goes to WD heterozygotes. These have genotype 7.
    # WD_pos = genotype_position_unordered(4,1);    [=7]
    females[7] += ceil(Int, 2*popsize[1,1]*xD/(1-2*xD));
    males[7] += ceil(Int, 2*popsize[2,1]*yD/(1-2*yD));
    # Now proportion of drive alleles satisfy x[4,1]=xD, y[4,1]=yD (approximately, as value is rounded).

    for t = 1:T-1
        mated_females = mating(females,males,pmat);
        # Migration does not occur in this algorithm
        # So skip to reproduction
        # Now for females and males in next generation
        totfem = sum(females);
        females = zeros(Int,n);
        males = zeros(Int,n);
        for kk = 1:n^2
            (female_genotype_offspring,male_genotype_offspring) = offspring_production(kk,mated_females[kk],totfem,Rm,K,fitness_vector,gamete_array);
            females += female_genotype_offspring;
            males += male_genotype_offspring;
        end

        popsize[:,t+1] = [sum(females);sum(males);sum(females+males)];
        # Each allele contributes 1 towards xnum. So an RR individual contributes 2 to xnum[2,t].
        # This means that sum(xnum[:,t]) = 2*sum(females). We will halve xnum, ynum at the end to correct for this.
        for i = 1:n
            xnum[genotype_list[i,1],t+1] += females[i];
            xnum[genotype_list[i,2],t+1] += females[i];
            ynum[genotype_list[i,1],t+1] = males[i];
            ynum[genotype_list[i,2],t+1] = males[i];
        end
        x[:,t+1] = xnum[:,t+1];
        y[:,t+1] = ynum[:,t+1];

        if maximum(x[:,t+1]) > 0
            x[:,t+1] = x[:,t+1]/sum(x[:,t+1]);
        end
        if maximum(y[:,t+1]) > 0
            y[:,t+1] = y[:,t+1]/sum(y[:,t+1]);
        end

        if x[2,t+1]+x[3,t+1] > test_threshold && y[2,t+1]+y[3,t+1] >test_threshold
            resistance = 1;
            τ = t;
        end
    end

    # Correct xnum, ynum as discussed above
    xnum = xnum/2;
    ynum = ynum/2;

    return (x,y,popsize,xnum,ynum,resistance,τ)
end

function fixedmig_TSR_mated(s,hD,hN,σ,Rm,K,xD,yD,ϵ,ν,μ,β,ξ,width::Int,height::Int,pmig,pmat,preexisting,T,test_threshold)
    # m=1 gRNA only for this program
    # No Gauss-Poisson hybrid for this
    Gauss = 0;
    #Set the number of alleles, r, and the number of genotypes, n = (1/2)r(r+1)
    r = 4;
    n = 10;
    # genotype_list[i,:] contains the alleles within genotype i.
    genotype_list = allele_comb_unordered(r);

    # First calculate gamete frequencies for the genotypes (see gamete_freq_array function for more details) and females fitnesses (see female_fitness funtion for details)
    gamete_array = gamete_freq_array(ϵ,ν,μ,β,ξ);
    fitness_vector = female_fitness(s,hD,hN,σ);

    # Calculate local carrying capacity
    K_perdeme = div(K,width*height,RoundUp);

    # At first all individuals have genotype 1 (homozygous W)
    females = zeros(Int,width,height,n);
    males = zeros(Int,width,height,n);
    for i = 1:width, j=1:height
        females[i,j,1] = rand(Poisson(0.5*K_perdeme));
        males[i,j,1] = rand(Poisson(0.5*K_perdeme));
    end

    mated_females = zeros(width,height,n^2); # Throw away unmated females to save memory

    # Now burn-in if wanted
    if preexisting == 1
        for t = 1:ceil(Int,5/σ)   # Longer burn-in time than initial program as we don't start from an estimate of the standing variation
            for i=1:width, j=1:height
                mated_females[i,j,:] = mating(females[i,j,:],males[i,j,:],pmat)[1:end-1];
            end
            ### Throw away females which did not find a mate.
            ###mated_females = mated_females[1:end-1];

            # Mated females migrate
            for ii = 1:n^2
                mated_females[:,:,ii] = square_deme_migration(mated_females[:,:,ii],pmig,Gauss);
            end
            # Now for females and males in next generation
            for i = 1:width, j = 1:height
                totfem = sum(females[i,j,:]);
                females[i,j,:] = zeros(Int,n);
                males[i,j,:] = zeros(Int,n);
                for kk = 1:n^2
                    (female_genotype_offspring,male_genotype_offspring) = offspring_production(kk,mated_females[i,j,kk],totfem,Rm,K_perdeme,fitness_vector,gamete_array);
                    females[i,j,:] += female_genotype_offspring;
                    males[i,j,:] += male_genotype_offspring;
                end
            end
        end
    end

    # Main phase of simulation
    # Will output (x,y,popsize,xnum,ynum,resistance,τ)
    xnum = zeros(width,height,r,T);
    ynum = zeros(width,height,r,T);
    x = zeros(width,height,r,T);
    y = zeros(width,height,r,T);
    popsize = zeros(3,width,height,T);
    # Initialise resistance, τ
    resistance = 0;
    τ = NaN;

    for i = 1:width, j = 1:height
        popsize[:,i,j,1] = [sum(females[i,j,:]);sum(males[i,j,:]);sum(females[i,j,:]+males[i,j,:])];
    end

    # Each allele contributes 1 towards xnum. So an RR individual contributes 2 to xnum[2,t].
    # This means that sum(xnum[:,t]) = 2*sum(females). We will halve xnum, ynum at the end to correct for this.
    for i = 1:width, j = 1:height, kk = 1:n
        xnum[i,j,genotype_list[kk,1],1] += females[i,j,kk];
        xnum[i,j,genotype_list[kk,2],1] += females[i,j,kk];
        ynum[i,j,genotype_list[kk,1],1] = males[i,j,kk];
        ynum[i,j,genotype_list[kk,2],1] = males[i,j,kk];
    end

    for i = 1:width, j = 1:height
        x[i,j,:,1] = xnum[i,j,:,1];
        y[i,j,:,1] = ynum[i,j,:,1];

        if maximum(x[i,j,:,1]) > 0
            x[i,j,:,1] = x[i,j,:,1]/sum(x[i,j,:,1]);
        end
        if maximum(y[i,j,:,1]) > 0
            y[i,j,:,1] = y[i,j,:,1]/sum(y[i,j,:,1]);
        end

        # Now add drive. For now, all drive goes to WD heterozygotes. These have genotype 7.
        # WD_pos = genotype_position_unordered(4,1);    [=7]
        females[i,j,7] += ceil(Int, 2*popsize[1,i,j,1]*xD/(1-2*xD));
        males[i,j,7] += ceil(Int, 2*popsize[2,i,j,1]*yD/(1-2*yD));
        # Now proportion of drive alleles satisfy x[4,1]=xD, y[4,1]=yD (approximately, as value is rounded).
    end

    for t = 1:T-1
        for i=1:width, j=1:height
            mated_females[i,j,:] = mating(females[i,j,:],males[i,j,:],pmat)[1:end-1];
        end
        ### Throw away females which did not find a mate.
        ###mated_females = mated_females[1:end-1];

        # Mated females migrate
        for ii = 1:n^2
            mated_females[:,:,ii] = square_deme_migration(mated_females[:,:,ii],pmig,Gauss);
        end
        # Now for females and males in next generation
        for i = 1:width, j = 1:height
            totfem = sum(females[i,j,:]);
            females[i,j,:] = zeros(Int,n);
            males[i,j,:] = zeros(Int,n);
            for kk = 1:n^2
                (female_genotype_offspring,male_genotype_offspring) = offspring_production(kk,mated_females[i,j,kk],totfem,Rm,K_perdeme,fitness_vector,gamete_array);
                females[i,j,:] += female_genotype_offspring;
                males[i,j,:] += male_genotype_offspring;
            end

            popsize[:,i,j,t+1] = [sum(females[i,j,:]);sum(males[i,j,:]);sum(females[i,j,:]+males[i,j,:])];
        end

        # Each allele contributes 1 towards xnum. So an RR individual contributes 2 to xnum[2,t].
        # This means that sum(xnum[:,t]) = 2*sum(females). We will halve xnum, ynum at the end to correct for this.
        for i = 1:width, j = 1:height, kk = 1:n
            xnum[i,j,genotype_list[kk,1],t+1] += females[i,j,kk];
            xnum[i,j,genotype_list[kk,2],t+1] += females[i,j,kk];
            ynum[i,j,genotype_list[kk,1],t+1] = males[i,j,kk];
            ynum[i,j,genotype_list[kk,2],t+1] = males[i,j,kk];
        end

        for i = 1:width, j = 1:height
            x[i,j,:,t+1] = xnum[i,j,:,t+1];
            y[i,j,:,t+1] = ynum[i,j,:,t+1];

            if maximum(x[i,j,:,t+1]) > 0
                x[i,j,:,t+1] = x[i,j,:,t+1]/sum(x[i,j,:,t+1]);
            end
            if maximum(y[i,j,:,t+1]) > 0
                y[i,j,:,t+1] = y[i,j,:,t+1]/sum(y[i,j,:,t+1]);
            end
        end

        if sum(xnum[:,:,2,t+1]) > test_threshold && sum(ynum[:,:,2,t+1]) >test_threshold
            resistance = 1;
            τ = t;
        end
    end

    return (x,y,popsize,xnum,ynum,resistance,τ)

end
