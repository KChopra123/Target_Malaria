function allele_comb(n)
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

function genotype_position(k1,k2)
    if k1 < k2
        (k1,k2) = (k2,k1);
    end
    return trunc(Int, 0.5*k1*(k1-1) + k2)
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
    # We count pairs of genotypes in the following way, as given by the function allele_comb above.
    # Consider n genotypes called 1,2,3,...,n. There are (1/2)n(n+1) unordered pairs.
    # The ordering of the set of pairs used is (1,1), (2,1), (2,2), (3,1), (3,2), (3,3), (4,1), (4,2), (4,3), (4,4), (5,1), ...
    # Including females which do not find a mate, there are (1/2)n(n+1)+1 entries.
    mated_females = zeros(Int, trunc(Int,0.5*n*(n+1)) + 1);
    #println("size of mated females vector is $(size(mated_females))")

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
    # This will gve a mated female with genotype (unordered) pair (i,j).
    # Unmated females in the entry mated_females[n+1] will obviously not reproduce.
    for i = 1:n
        v = rand(Multinomial(females[i],prob_vec));
        for j = 1:n
            mated_females[genotype_position(i,j)] += v[j];
        end
        mated_females[end] += v[end];
    end

    return mated_females

end