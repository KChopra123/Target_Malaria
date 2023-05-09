function square_deme_size(D,pmig,area_approx)
# Approximates the size of deme needed to *roughly* simulate a Brownian motion
# with diffusivity D and probability of migration pmig

    # Can analytically approximate by setting numerical = 0
    if area_approx == 1
        r = 2*sqrt(-D * log(pmig));
        x = sqrt(pi) * r;
    # Or can use the in-built error function to find the size more exactly
    else
        x = sqrt(D) * 4 * erfinv(sqrt(1 - pmig));
    end

    # Output the deme size length, x
    return x
end
