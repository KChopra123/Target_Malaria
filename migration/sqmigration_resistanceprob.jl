function sqmigration_resistanceprob(L,m,s,h,hN,σ,Rm,K,xD,yD,ϵ,ν,μ,β,ξ,preexisting,T,test_threshold,width,height,Gauss,D,pmig,area_approx)
    # Runs L simulations of sqmigration_TSR with the above parameters.
    # The burn-in is quite computationally expensive, so will be done once
    # and used for every simulation.
    #Set proportion of males
    ψ = 0.5;

    N = round(Int128,K);
    if Gauss == 0 && log2(N) > 63
        print("Error - if population > 2^63 then must use Gauss-Poisson Approximation")
        return -2
    end

    x = square_deme_size(D,pmig,area_approx);

    # To fit the area to square demes, we use the smallest number of demes in
    # each direction needed to cover the space
    acrm = div(width,x,RoundUp);
    acrm = round(Int,acrm);
    upn = div(height,x,RoundUp);
    upn = round(Int,upn);
    println("checkpoint 1")
    A = zeros(Int128,(acrm,upn));

    # Now we distribute the individuals amongst the demes.
    # For now, assume uniform distribution. This can be changed later.
    K_perdeme = div(K,acrm*upn,RoundUp);
    for i = 1:acrm
        for j = 1:upn
            A[i,j] = K_perdeme;
        end
    end
    println("checkpoint 2")

    # Sort population per deme into males and females
    Nm = zeros(Int128,(acrm,upn));
    Nf = zeros(Int128,(acrm,upn));

    if Gauss == 0
        for i = 1:acrm
            for j = 1:upn
                Nm[i,j] = rand(Binomial(A[i,j],ψ));
                Nf[i.j] = A[i,j] - Nm[i,j];
            end
        end
    else
        for i = 1:acrm
            for j = 1:upn
                vvvv = GaussPoissonHybrid_mnrnd(A[i,j],[ψ;1-ψ]);
                Nm[i,j] = round(Int128(vvvv[1]));
                Nf[i,j] = A[i,j] - Nm[i,j];
            end
        end
    end

    N = Nf + Nm;
    println("checkpoint 3")

    # Density dependent parameter from Beverton-Holt model of population
    # dynamics
    # Note that we use the carrying capacity of *each deme*
    α = K_perdeme/(Rm - 1);

    # Dominance coefficients
    hD = h;

    if length(σ) == 1
        σ = σ*fill(1.0,(m,1));
    elseif length(σ) != m
        print("Error - σ should be a scalar of vector of length m")
        return -3
    end
    println("checkpoint 4")

    # Set number of alleles (this would be difficult to change anyway in
    # practice)
    n = 4;

    # calculate number of (ordered) haplotypes
    na  = (n-1)^m + 1;

    #calculate fitness matrices
    Wm = fill(1.0,(na,na));
    Wf = femalefitnessmatrix(m,n,s,hD,hN,σ);

    #Calculate mutation probaility matrix. M[i,j] is the probability of de novo
    #mutation haplotype j -> haplotype i
    mutmatrix = mutationmatrix(m,n,μ,ξ);

    #Calculate drive conversion matrix. Suppose an individual has genotype
    #[j,drive]. Then K[i,j]+δ_{ij} is the probability of the conversion j->i via
    #drive and NHEJ mutations. So K[i,j] is the mean relative change in
    #frequency j->i.
    κ = drivematrix(m,n,ϵ,ν,β);
    println("checkpoint 5")

    ## Initial frequencies
    x0 = zeros((acrm,upn,na));
    y0 = zeros((acrm,upn,na));
    for i = 1:acrm
        for j = 1:upn
            x0[i,j,1] = 1.0;
            y0[i,j,1] = 1.0;
        end
    end
    println("checkpoint 6")

    ## If preexisting
    if preexisting == 1
        #Calculate burn in time
        if minimum(σ) <= 0
            print("Error - for preexisting = 1, resistance alleles must always incur fitness cost")
            return
        end

        Tburn = round(1/minimum(σ));

        # Now initialise as in original code
        xxx0 = zeros(na);
        yyy0 = zeros(na);
        xxx0[1] = 1.0;
        yyy0[1] = 1.0;
        println("checkpoint 7")

        indm = mutmatrix.>0;
        for i = 1:na
            for j = 2:na
                indm[i,j] = 0;
            end
        end
        indm[1,1] = 0;

        muin = mutmatrix[indm];
        #print("muin = ")
        #println(muin)

        ffff = diagm((Wf[indm].-1)./2);
        #print("ffff = ")
        #println(ffff)
        #println(ffff)

        indm2 = mutmatrix.>0;
        indm2[1,1] = 0;
        for i = 2:na
            for j = 2:na
                if indm2[i,1] > 0 && indm2[j,1] > 0
                    indm2[i,j] = 1;
                else
                    indm2[i,j] = 0;
                end
            end
        end

        for i = 1:na
            indm2[i,1] = 0;
            indm2[1,i] = 0;
        end

        modmut1 = mutmatrix[indm2];
        lengthvec = length(modmut1);
        if lengthvec != length(ffff);
            println("Error with preexisting")
            println(lengthvec)
            println(length(ffff))
        end
        lengthvec = round(Int,sqrt(lengthvec));
        modmut2 = zeros((lengthvec,lengthvec));
        for i = 1:lengthvec
            for j = 1:lengthvec
                modmut2[i,j] = modmut1[lengthvec*(j-1)+i];
                if i == j
                    modmut2[i,j] += -1;
                end
            end
        end
        #print("modmut2 = ")
        #println(modmut2)

        indmvec = BitArray(undef,na);
        for i = 1:na
            indmvec[i] = indm[i,1];
        end

        xxx0[indmvec] = -inv(modmut2+ffff)*muin;
        yyy0[indmvec] = -inv(modmut2+ffff)*muin;
        #println(x0)
        #println(y0)

        xxx0[na] = 0;
        yyy0[na] = 0;
        xxx0[1] = 1 - sum(xxx0[2:end]);
        yyy0[1] = 1 - sum(yyy0[2:end]);
        println("checkpoint 8")
        print("xxx0 = ")
        println(xxx0)

        for i = 1:acrm
            for j = 1:upn
                x0[i,j,:] = xxx0;
                y0[i,j,:] = yyy0;
            end
        end
        println("checkpoint 9, burn-in commencing")

        for tt = 1:Tburn
            println(tt)
            (x0,y0,Nf,Nm,N) = sqmignewfrequency(Rm,α,x0,y0,μ,Wf,Wm,na,ψ,mutmatrix,κ,Nf,Nm,Gauss,pmig);
        end
        println("checkpoint 10, burn-in complete")
    end

    # Now add drive
    for i = 1:acrm
        for j = 1:upn
            x0[i,j,na] = xD;
            x0[i,j,1] = x0[i,j,1] - xD;
            y0[i,j,na] = yD;
            y0[i,j,1] = y0[i,j,1] - yD;

            if x0[i,j,1] < 0 || y0[i,j,1] < 0
                println("Error - not enough wild alleles to add drive in deme")
                println([i,j])
                return 5
            end
        end
    end

    # The initial frequencies are now x0 and y0 for every simulation
    xcurrent = x0;
    ycurrent = y0;

    x = zeros(L,acrm,upn,na,T);
    y = zeros(L,acrm,upn,na,T);
    time = zeros(Int,T);

    popsize = zeros(Int128,(L,3,acrm,upn,T));

    testallele = 2;
    resistance = zeros(Int,L);
    τ = zeros(Int,L);
    resistancepositions = resistancehaplotypes(m,n,[2,3]);

    Nf0 = Nf;
    Nm0 = Nm;

    for l = 1:L
        print("simulation number ")
        println(l)
        x[l,:,:,:,1] = x0;
        y[l,:,:,:,1] = y0;
        popsize[l,1,:,:,1] = Nf0;
        popsize[l,2,:,:,1] = Nm0;
        popsize[l,3,:,:,1] = Nf0+Nm0;

        for t = 1:T-1
            println(t)
            time[t+1] = t;
            newfrequencies = sqmignewfrequency(Rm,α,xcurrent,ycurrent,μ,Wf,Wm,na,ψ,mutmatrix,κ,Nf,Nm,Gauss,pmig);
            if newfrequencies[1] == 55 || newfrequencies[1] == 56
                (~,iii,jjj,kkkk) = newfrequencies;
                println(x[l,iii,jjj,:,1:t+1])
                println(y[l,iii,jjj,:,1:t+1])
                println(popsize[l,:,iii,jjj,1:t+1])
                return (x[l,iii,jjj,:,1:t+1],y[l,iii,jjj,:,1:t+1],popsize[l,:,iii,jjj,1:t+1])
            else
                (xcurrent,ycurrent,Nf,Nm,N) = newfrequencies;
            end

            x[l,:,:,:,t+1] = xcurrent;
            y[l,:,:,:,t+1] = ycurrent;
            popsize[l,1,:,:,t+1] = Nf;
            popsize[l,2,:,:,t+1] = Nm;
            popsize[l,3,:,:,t+1] = Nf+Nm;

            # Test for resistance
            # Discuss exactly how do do this with supervisor. For now, do the following
            resistfreqf = zeros(acrm,upn);
            resistfreqm = zeros(acrm,upn);
            for i = 1:acrm
                for j = 1:upn
                    resistfreqf[i,j] = dot(xcurrent[l,i,j,:],resistancepositions);
                    resistfreqm[i,j] = dot(ycurrent[l,i,j,:],resistancepositions);
                end
            end
            scalarresistf = mean(resistfreqf);
            scalarresistm = mean(resistfreqm);
            if scalarresistf >= test_threshold  && scalarresistm >= test_threshold
                if resistance[l] == 0
                    τ[l] = t;
                end
                resistance[l] = 1;
            end
        end

        if resistance[l] == 0
            τ[l] = NaN;
        end

    end

    xnum = zeros(size(x));
    ynum = zeros(size(y));
    for i = 1:na
        for t = 1:T
            xnum[:,:,i,t] = x[:,:,i,t].*popsize[1,:,:,t];
            ynum[:,:,i,t] = y[:,:,i,t].*popsize[2,:,:,t];
        end
    end

    return (x,y,popsize,xnum,ynum,resistance,τ)

end
