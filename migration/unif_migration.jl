function unif_migration(A)
    # A is a vector of population sizes
    # Let n denote the length of A. So A = (A_i) for 1≦i≦n.
    # Each individual is uniformly randomly placed in one of the n demes.

    n = length(A);
    S = sum(A);

    A = rand(Multinomial(S,n));

    return A
end
