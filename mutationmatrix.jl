function mutationmatrix(m,n,μ,ξ)
    #Outputs mutmatrix. mutmatrix[i,j] denotes the probability of de novo mutation from
    #haplotype j to haplotype i.

    #First calculate the number of (ordered) haplotypes
    na = (n-1)^m + 1;

    #muvec = [1 - muSite; xi * muSite; (1 - xi) * muSite];

    #uR = [0; 1; 0];
    #uN = [0; 0; 1];

    #singlesitemut is the matrix of mutation probabilities j->i for m=1
    singlesitemut = [1-μ 0.0 0.0 0.0;ξ*μ 1.0 0.0 0.0;(1-ξ)*μ 0.0 1.0 0.0;0.0 0.0 0.0 1.0];

    #Probability of mutation from j to i
    mutmatrix = zeros(na,na);

    indexarray = linindices(m,n);

    for i = 1:na-1
        for j = 1:na-1
            for kk = 1:m
                ii = indexarray[i,kk];
                jj = indexarray[j,kk];
                #all entries are initially zero, so have to set first value
                if kk == 1
                    mutmatrix[i,j] = singlesitemut[ii,jj];
                else
                    mutmatrix[i,j] = mutmatrix[i,j]*singlesitemut[ii,jj];
                end
            end
        end
    end

    mutmatrix[na,na] = 1.0;

    #subtract identity matrix??

    return mutmatrix;

end
