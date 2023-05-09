function TSR_resistance_threshold(s,h,hN,σ,n_r,Rm,K,xD,yD,ϵ,ν,width,height,pmig,T,test_threshold,Gauss)
    # Adds single (female) functional resistant allele to a uniformly random deme
    m=1;
    μ = 0;
    β = 0;
    ξ = 0;
    terminate_resistance = 1;
    #Set proportion of males
    ψ = 0.5;

    N = round(Int128,K);
    if Gauss == 0 && log2(N) > 63
        print("Error - if population > 2^63 then must use Gauss-Poisson Approximation")
        return -2
    end

    # To fit the area to square demes, we use the smallest number of demes in
    # each direction needed to cover the space
    acrm = width;
    acrm = round(Int,acrm);
    upn = height;
    upn = round(Int,upn);
    #println("checkpoint 1")
    A = zeros(Int128,(acrm,upn));

    # Now we distribute the individuals amongst the demes.
    # For now, assume uniform distribution. This can be changed later.
    K_perdeme = div(K,acrm*upn,RoundUp);
    for i = 1:acrm
        for j = 1:upn
            A[i,j] = K_perdeme;
        end
    end
    #println("checkpoint 2")

    # Sort population per deme into males and females
    Nm = zeros(Int128,(acrm,upn));
    Nf = zeros(Int128,(acrm,upn));

    if Gauss == 0
        for i = 1:acrm
            for j = 1:upn
                Nm[i,j] = rand(Binomial(A[i,j],ψ));
                Nf[i,j] = A[i,j] - Nm[i,j];
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
    #println("checkpoint 3")

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
    #println("checkpoint 4")

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
    #println("checkpoint 5")

    ## Initial frequencies
    x0 = zeros((acrm,upn,na));
    y0 = zeros((acrm,upn,na));
    for i = 1:acrm
        for j = 1:upn
            x0[i,j,1] = 1.0;
            y0[i,j,1] = 1.0;
        end
    end

    # Add in resistance
    for i = 1:n_r
        x_coord = rand(1:acrm);
        y_coord = rand(1:upn);
        r_0 = 1/Nf[x_coord,y_coord];
        x0[x_coord,y_coord,2] = x0[x_coord,y_coord,2] + r_0;
        x0[x_coord,y_coord,1] = 1 - x0[x_coord,y_coord,2];
    end

    #println("checkpoint 6")

    for i = 1:acrm
        for j = 1:upn
            if i > 1 || j > 1
                x0[i,j,na] = xD;
                x0[i,j,1] = x0[i,j,1] - xD;
                y0[i,j,na] = yD;
                y0[i,j,1] = y0[i,j,1] - yD;
            end

            if x0[i,j,1] < 0 || y0[i,j,1] < 0
                println("Error - not enough wild alleles to add drive in deme")
                println([i,j])
                return 5
            end
        end
    end

    xcurrent = x0;
    ycurrent = y0;
    #println("checkpoint 13")

    # x,y will keep track over time of fe/male haplotype frequencies,
    # respectively
    #x = zeros(acrm,upn,na,T);
    #y = zeros(acrm,upn,na,T);
    #time = zeros(Int,T);
    #x[:,:,:,1] = xcurrent;
    #y[:,:,:,1] = ycurrent;
    #println("checkpoint 14")

    # popsize keeps track of female, male and total populations
    #popsize = zeros(Int128,(3,acrm,upn,T));
    #popsize[1,:,:,1] = Nf;
    #popsize[2,:,:,1] = Nm;
    #popsize[3,:,:,1] = Nf+Nm;
    #println("checkpoint 15")

    # We want to keep track of the resistant allele
    testallele = 2;
    resistance = 0;
    resistancepositions = resistancehaplotypes(m,n,[2]);
    #println("checkpoint 16")

    for t = 1:T-1
        print("$t ")
        #time[t+1] = t;
        #(xcurrent,ycurrent,Nf,Nm,N) = sqmignewfrequency(Rm,α,xcurrent,ycurrent,μ,Wf,Wm,na,ψ,mutmatrix,κ,Nf,Nm,Gauss,pmig);
        newfrequencies = sqmignewfrequency(Rm,α,xcurrent,ycurrent,μ,Wf,Wm,na,ψ,mutmatrix,κ,Nf,Nm,Gauss,pmig);
        if newfrequencies[1] == 55 || newfrequencies[1] == 56
            (~,iii,jjj,kkkk) = newfrequencies;
            #println(x[iii,jjj,:,1:t+1])
            #println(y[iii,jjj,:,1:t+1])
            #println(popsize[:,iii,jjj,1:t+1])
            return (x[iii,jjj,:,1:t+1],y[iii,jjj,:,1:t+1],popsize[:,iii,jjj,1:t+1],resistance,τ)
        else
            (xcurrent,ycurrent,Nf,Nm,N) = newfrequencies;
        end

        #x[:,:,:,t+1] = xcurrent;
        #y[:,:,:,t+1] = ycurrent;
        #popsize[1,:,:,t+1] = Nf;
        #popsize[2,:,:,t+1] = Nm;
        #popsize[3,:,:,t+1] = Nf+Nm;

        # Test for resistance
        # Discuss exactly how do do this with supervisor. For now, do the following
        resistfreqf = zeros(acrm,upn);
        resistfreqm = zeros(acrm,upn);
        for i = 1:acrm
            for j = 1:upn
                resistfreqf[i,j] = dot(xcurrent[i,j,:],resistancepositions);
                resistfreqm[i,j] = dot(ycurrent[i,j,:],resistancepositions);
            end
        end
        #scalarresistf = mean(resistfreqf);
        #scalarresistm = mean(resistfreqm);
        scalarresistf = 0;
        scalarresistm = 0;
        for i = 1:acrm
            for j = 1:upn
                scalarresistf += resistfreqf[i,j]*Nf[i,j];
                scalarresistm += resistfreqm[i,j]*Nm[i,j];
            end
        end
        scalarresistf = scalarresistf;
        scalarresistm = scalarresistm;
        if scalarresistf == 0 && scalarresistm == 0
            #resistance = -1;
            #τ = t;
            return (resistance,NaN)
        end
        if scalarresistf >= test_threshold  && scalarresistm >= test_threshold
            if resistance == 0
                τ = t;
            end
            resistance = 1;
            if terminate_resistance == 1
                #write_matfile("long_sim.mat"; x=x, y=y, popsize=popsize, resistance=resistance)
                return (resistance,τ)
            end
        end

    end
    if resistance == 0;
        τ = NaN;
    end
    return (resistance,τ)
end
