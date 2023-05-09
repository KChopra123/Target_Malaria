function square_deme_migration(A,pmig,Gauss)
# The entry A[i,j] denotes the initial population of individuals in deme [i,j].
# Individuals migrate between dquare demes
# Each individual migrates with probability pmig. The destination is chosen
# uniformly randomly from the neighbouring demes
    acrm = size(A)[1];
    upn = size(A)[2];

    # This program only works for at least 2x2 matrices
    if acrm<2 || upn<2
        print("Error - need at least 2x2 grid for square demes")
        return 1
    end
    #println("α")

    # Create arrays of the number of individuals emigrating from each deme.
    # For example, remain[i,j] denotes the number of individuals moving
    # [i,j]->[i,j], whereas down[i,j]  represents the number moving
    # [i,j]->[i,j-1].
    if Gauss == 0
        remain = zeros(Int,(acrm,upn));
        up = zeros(Int,(acrm,upn));
        down = zeros(Int,(acrm,upn));
        left = zeros(Int,(acrm,upn));
        right = zeros(Int,(acrm,upn));
    else
        remain = zeros(Int128,(acrm,upn));
        up = zeros(Int128,(acrm,upn));
        down = zeros(Int128,(acrm,upn));
        left = zeros(Int128,(acrm,upn));
        right = zeros(Int128,(acrm,upn));
    end
    #println("β")

    for i = 1:acrm
        for j = 1:upn
            if Gauss == 0
                # Need population to have type Int
                A_temp = zeros(Int,(acrm,upn));
                for i = 1:acrm
                    for j = 1:upn
                        A_temp[i,j] = trunc(Int,A[i,j]);
                    end
                end
                A = A_temp;
                # First consider interior demes
                if 1<i<acrm && 1<j<upn
                    #println("γ")
                    v = rand(Multinomial(A[i,j],[1-pmig;pmig/4;pmig/4;pmig/4;pmig/4]));
                    #println("δ")
                    remain[i,j] = v[1];
                    right[i,j] = v[2];
                    up[i,j] = v[3];
                    left[i,j] = v[4];
                    down[i,j] = v[5];

                # Now non-corner demes on the left edge
                elseif i == 1 && 1<j<upn
                    v = rand(Multinomial(A[i,j],[1-(3/4)*pmig;pmig/4;pmig/4;pmig/4]));

                    remain[i,j] = v[1];
                    right[i,j] = v[2];
                    up[i,j] = v[3];
                    # no left movement
                    down[i,j] = v[4];

                # Now right edge
                elseif i == acrm && 1<j<upn
                    v = rand(Multinomial(A[i,j],[1-(3/4)*pmig;pmig/4;pmig/4;pmig/4]));

                    remain[i,j] = v[1];
                    up[i,j] = v[2];
                    left[i,j] = v[3];
                    down[i,j] = v[4];

                # Now bottom edge
                elseif 1<i<acrm && j == 1
                    v = rand(Multinomial(A[i,j],[1-(3/4)*pmig;pmig/4;pmig/4;pmig/4]));

                    remain[i,j] = v[1];
                    right[i,j] = v[2];
                    up[i,j] = v[3];
                    left[i,j] = v[4];

                elseif 1<i<acrm && j == upn
                    v = rand(Multinomial(A[i,j],[1-(3/4)*pmig;pmig/4;pmig/4;pmig/4]));

                    remain[i,j] = v[1];
                    right[i,j] = v[2];
                    left[i,j] = v[3];
                    down[i,j] = v[4];

                elseif i == 1 && j == 1
                    #println("γ")
                    v = rand(Multinomial(A[i,j],[1-(1/2)*pmig;pmig/4;pmig/4]));
                    #println("δ")
                    remain[i,j] = v[1];
                    right[i,j] = v[2];
                    up[i,j] = v[3];

                elseif i == acrm && j == 1
                    v = rand(Multinomial(A[i,j],[1-(1/2)*pmig;pmig/4;pmig/4]));
                    remain[i,j] = v[1];
                    up[i,j] = v[2];
                    left[i,j] = v[3];

                elseif i == 1 && j == upn
                    v = rand(Multinomial(A[i,j],[1-(1/2)*pmig;pmig/4;pmig/4]));
                    remain[i,j] = v[1];
                    right[i,j] = v[2];
                    down[i,j] = v[3];

                elseif i == acrm && j == upn
                    v = rand(Multinomial(A[i,j],[1-(1/2)pmig;pmig/4;pmig/4]));
                    remain[i,j] = v[1];
                    left[i,j] = v[2];
                    down[i,j] = v[3];
                else
                    print("Error calculating emigrant numbers")
                    return 2
                end
            else
                # Need population to have type Int128
                A_temp = zeros(Int128,(acrm,upn));
                for i = 1:acrm
                    for j = 1:upn
                        A_temp[i,j] = trunc(Int,A[i,j]);
                    end
                end
                A = A_temp;
                # First consider interior demes
                if 1<i<acrm && 1<j<upn
                    #println("γ")
                    v = GaussPoissonHybrid_mnrnd(A[i,j],[1-pmig;pmig/4;pmig/4;pmig/4;pmig/4]);
                    #println("δ")
                    remain[i,j] = v[1];
                    right[i,j] = v[2];
                    up[i,j] = v[3];
                    left[i,j] = v[4];
                    down[i,j] = v[5];

                # Now non-corner demes on the left edge
                elseif i == 1 && 1<j<upn
                    v = GaussPoissonHybrid_mnrnd(A[i,j],[1-pmig;pmig/3;pmig/3;pmig/3]);

                    remain[i,j] = v[1];
                    right[i,j] = v[2];
                    up[i,j] = v[3];
                    # no left movement
                    down[i,j] = v[4];

                # Now right edge
                elseif i == acrm && 1<j<upn
                    v = GaussPoissonHybrid_mnrnd(A[i,j],[1-pmig;pmig/3;pmig/3;pmig/3]);

                    remain[i,j] = v[1];
                    up[i,j] = v[2];
                    left[i,j] = v[3];
                    down[i,j] = v[4];

                # Now bottom edge
                elseif 1<i<acrm && j == 1
                    v = GaussPoissonHybrid_mnrnd(A[i,j],[1-pmig;pmig/3;pmig/3;pmig/3]);

                    remain[i,j] = v[1];
                    right[i,j] = v[2];
                    up[i,j] = v[3];
                    left[i,j] = v[4];

                elseif 1<i<acrm && j == upn
                    v = GaussPoissonHybrid_mnrnd(A[i,j],[1-pmig;pmig/3;pmig/3;pmig/3]);

                    remain[i,j] = v[1];
                    right[i,j] = v[2];
                    left[i,j] = v[3];
                    down[i,j] = v[4];

                elseif i == 1 && j == 1
                    #println("γ")
                    v = GaussPoissonHybrid_mnrnd(A[i,j],[1-pmig;pmig/2;pmig/2]);
                    #println(v)
                    #println("δ")
                    remain[i,j] = v[1];
                    #println("δ.5")
                    right[i,j] = v[2];
                    up[i,j] = v[3];
                    #println("ϵ")

                elseif i == acrm && j == 1
                    v = GaussPoissonHybrid_mnrnd(A[i,j],[1-pmig;pmig/2;pmig/2]);
                    remain[i,j] = v[1];
                    up[i,j] = v[2];
                    left[i,j] = v[3];

                elseif i == 1 && j == upn
                    v = GaussPoissonHybrid_mnrnd(A[i,j],[1-pmig;pmig/2;pmig/2]);
                    remain[i,j] = v[1];
                    right[i,j] = v[2];
                    down[i,j] = v[3];

                elseif i == acrm && j == upn
                    v = GaussPoissonHybrid_mnrnd(A[i,j],[1-pmig;pmig/2;pmig/2]);
                    remain[i,j] = v[1];
                    left[i,j] = v[2];
                    down[i,j] = v[3];
                else
                    print("Error calculating emigrant numbers")
                    return 2
                end
            end
        end
    end

    for i = 1:acrm
        if down[i,1] != 0
            print("Error calculating emigrant numbers")
            return 3
        elseif up[i,upn] != 0
            print("Error calculating emigrant numbers")
            return 4
        end
    end

    for j = 1:upn
        if left[1,j] != 0
            print("Error calculating emigrant numbers")
            return 5
        elseif right[acrm,j] != 0
            print("Error calculating emigrant numbers")
            return 6
        end
    end

    # Now calculate population distribution after migration, B[i,j]
    if Gauss == 0
        B = zeros(Int,(acrm,upn));
    else
        B = zeros(Int128,(acrm,upn));
    end

    for i = 1:acrm
        for j = 1:upn
            if 1<i<acrm && 1<j<upn
                B[i,j] = remain[i,j] + right[i-1,j] + up[i,j-1] + left[i+1,j] + down[i,j+1];
            elseif i == 1 && 1<j<upn
                B[i,j] = remain[i,j] + up[i,j-1] + left[i+1,j] + down[i,j+1];
            elseif i == acrm && 1<j<upn
                B[i,j] = remain[i,j] + right[i-1,j] + up[i,j-1] + down[i,j+1];
            elseif 1<i<acrm && j == 1
                B[i,j] = remain[i,j] + right[i-1,j] + left[i+1,j] + down[i,j+1];
            elseif 1<i<acrm && j == upn
                B[i,j] = remain[i,j] + right[i-1,j] + up[i,j-1] + left[i+1,j];
            elseif i == 1 && j == 1
                B[i,j] = remain[i,j] + left[i+1,j] + down[i,j+1];
            elseif i == acrm && j == 1
                B[i,j] = remain[i,j] + right[i-1,j] + down[i,j+1];
            elseif i == 1 && j == upn
                B[i,j] = remain[i,j] + up[i,j-1] + left[i+1,j];
            elseif i == acrm && j == upn
                B[i,j] = remain[i,j] + right[i-1,j] + up[i,j-1];
            end
        end
    end

    return B

end
