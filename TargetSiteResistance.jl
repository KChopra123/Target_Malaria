function TargetSiteResistance(m,s,h,hN,σ,Rm,K,xD,yD,ϵ,ν,μ,β,ξ,preexisting,T,test_threshold,Gauss,Cswitch,plotfig,plotmale,vis)
# m = number of target sites
# s = homozygous fitness cost of drive
# s*h = fitness cost of heterozygous W/D individual
# s*hN = fitness cost of heterozygous W/N or R/D individual
# σ_i = fitness cost of R allele at target site i, 1≦i≦m
# Rm =
# K =
# xD,yD = initial proportion of D...D haplotypes in females/males respectively
# ϵ = efficiency of drive cleavage per target site
# ν = non-homologous end-joining (NHEJ) rate per generation per individual per
# target site
# μ = mutation rate per cleavage target site per generation per individual per
# target site (e.g. if a target site contains 10bp and bp mutation rate is 1e-9,
# μ=1e-8)
# β = fraction of functional resistant NHEJ mutations
# ξ = fraction of functional resistant single nucleotide mutations
# preexisting = 1,0 to run a burn in period of 1/σ and set initial allele
# frequencies or not, respectively
# T = length of simulations in generations

    #Set proportion of males
    ψ = 0.5;

    N=round(Int128,K);#Make this better later
    if Gauss == 0 && log2(N) > 63
        print("Error - if population > 2^63 then must use Gauss-Poisson Approximation")
        return 2
    end
    #Nm=ψ*N;
    #Nf=N-Nm;
    if Gauss == 0
        Nm = rand(Binomial(N,ψ));
        Nf = N - Nm;
    else
        vvv = GaussPoissonHybrid_mnrnd(N,[ψ;1-ψ]);
        Nm = round(Int128,vvv[1]);
        Nf = round(Int128,vvv[2]);
    end
    #Density dependent parameter from Beverton-Holt model of population dynamics
    α = K/(Rm - 1);

    #Dominance coefficients
    hD = h; #Why not just input hD instead of h and avoid this???

    if length(σ) == 1
        σ = σ*fill(1.0,(m,1));
    elseif length(σ) != m
        print("Error - σ should be scalar or vector of length m")
        return
    end

    #Set number of alleles
    n = 4;

    #Calculate number of (ordered) haplotypes
    na = (n-1)^m + 1;

    #Calculate fitness matrices
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

    ## Initial frequencies
    x0 = zeros(na);
    x0[1] = 1.0;
    y0 = zeros(na);
    y0[1] = 1.0;

    ##If preexisting
    if preexisting == 1 && μ > 0
        #Calculate burn in time
        if minimum(σ) <= 0
            print("Error - for preexisting = 1, resistance alleles must always incur fitness cost")
            return
        end

        Tburn = round(1/minimum(σ));

        # Now initialise as in original code
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

        if det(modmut2+ffff) != 0
        x0[indmvec] = -inv(modmut2+ffff)*muin;
        y0[indmvec] = -inv(modmut2+ffff)*muin;
        else
            Tburn = Tburn*10;
        end
        #println(x0)
        #println(y0)

        x0[na] = 0;
        y0[na] = 0;
        x0[1] = 1 - sum(x0[2:end]);
        y0[1] = 1 - sum(y0[2:end]);
        #print("x0 = ")
        #println(x0)

        for tt = 1:Tburn
            #m
            (x0,y0,Nf,Nm,N) = newfrequency(Rm,α,x0,y0,μ,Wf,Wm,na,ψ,mutmatrix,κ,Nf,Nm,Gauss);
            if (x0,y0,Nf,Nm,N) == (-1,-1,-1,-1,-1)
                println("error occured during burn in, time = ")
                print(tt)
                return
            end
        end
        #print(x0)
        #print(y0)
        #print(N)
    end

    ## Now run simulation, keeping track of frequencies and population size
    # First, add drive to initial frequencies
    x0[na] = xD;
    x0[1] = x0[1] - xD;
    y0[na] = yD;
    y0[1] = y0[1] - yD;

    if x0[1] < 0 || y0[1] < 0
        print("error - not enough wild alleles to add drive")
        return 5
    end

    xcurrent = x0;
    ycurrent = y0;
    # x,y will keep track over time of fe/male haplotype frequencies,
    # respectively
    x = zeros((na,T));
    y = zeros((na,T));
    time = zeros(Int,T);
    x[:,1] = xcurrent;
    y[:,1] = ycurrent;

    #popsize keeps track over time of female, male and total populations
    popsize = zeros(Int128, (3,T));
    popsize[:,1] = [Nf;Nm;N];

    #We want to keep track of the resistant allele
    testallele = 2;
    resistance = 0;
    resistancepositions = resistancehaplotypes(m,n,[2,3]);
    #Calculate the number of resistance alleles in each haplotype
    indexarray = linindices(m,n);

    for t = 1:T-1
        time[t+1] = t;
        (xcurrent,ycurrent,Nf,Nm,N) = newfrequency(Rm,α,xcurrent,ycurrent,μ,Wf,Wm,na,ψ,mutmatrix,κ,Nf,Nm,Gauss);
        if (xcurrent,ycurrent,Nf,Nm,N) == (-1,-1,-1,-1,-1)
            print("error occured during main phase, time = ")
            print(t)
            return(x,y,popsize)
        end
        x[:,t+1] = xcurrent;
        y[:,t+1] = ycurrent;
        popsize[:,t+1] = [Nf;Nm;N];

        # Test for resistance

        resistfreqf = dot(xcurrent,resistancepositions);
        resistfreqm = dot(ycurrent,resistancepositions);
        #resistancefreq = (1/N)*(Nf*resistfreqf + Nm*resistfreqm);
        if resistfreqf >= test_threshold && resistfreqm >= test_threshold
            if resistance == 0
                τ = t;
                #print(xcurrent)
                #print(ycurrent)
                #print(resistfreqf)
                #print(resistfreqm)
                #print(resistancepositions)
            end
            resistance = 1;
        end

    end

    if resistance == 0
        τ = NaN;
    end

    #recovery = 0;
    #if popsize[3,T] >= 0.1*K
    #    recovery = 1;
    #end

    return (x,y,popsize,resistance,τ)
end
