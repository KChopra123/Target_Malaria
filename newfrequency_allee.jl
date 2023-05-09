function newfrequency_allee(Rm,α,xold,yold,μ,Wf,Wm,na,ψ,mutmatrix,κ,Nf,Nm,pmat,Gauss)

    # Recover original carrying capacity. Define new constant.
    K = α * (Rm - 1);

    #if either sex has died out, there will be no reproduction
    if (mean(xold) == 0 || mean(yold) == 0) || (Nf == 0 || Nm == 0)
        xnew = zeros(na);
        ynew = zeros(na);
        return (xnew,ynew,0,0,0)
    end
    #First calculate average fitnesses
    meanf = transpose(xold)*Wf*yold;
    meanm = transpose(xold)*Wm*yold;
    if meanf == 0 || meanm == 0
        xnew = zeros(na);
        ynew = zeros(na);
        return (xnew,ynew,0,0,0)
    end

    #Now calculate the new frequencies obtained by random mating and ne novo
    #mutation. At this stage we NEGLECT DRIVE.
    xnewnodrive = mutmatrix * (((Wf*yold).*xold + (Wf*xold).*yold)/(2*meanf));
    ynewnodrive = mutmatrix * (((Wm*xold).*yold + (Wm*yold).*xold)/(2*meanm));

    #Now correct for drive
    xdriveterm = zeros(na);
    ydriveterm = zeros(na);

    for i = 1:na
        for j = 1:na
            #Check it's OK to allow all haplotypes. #Yes, I think it is.
            #H double counts drive homozygotes, but this is OK since κ zeros
            #them anyway
            #This has been modified from paper. Multiply by mutmatrix[j,j]
            #instead of (1-μ)
            H = xold[j]*yold[na] + xold[na]*yold[j];
            xdriveterm[i] += 0.5*Wf[j,na]*mutmatrix[j,j]*κ[i,j]*H/meanf;
            ydriveterm[i] += 0.5*Wm[j,na]*mutmatrix[j,j]κ[i,j]*H/meanm;
        end
        #print("xdriveterm = ")
        #println(xdriveterm)
        #print("ydriveterm = ")
        #println(ydriveterm)
    end

    xoverall = xnewnodrive + xdriveterm;
    if minimum(xoverall) < 0
        print("error - new expected female frequency has negative entry")
        return (-1,-1,-1,-1,-1)
    end
    yoverall = ynewnodrive + ydriveterm;
    if minimum(yoverall) < 0
        print("error - new expected male frequency has negative entry")
        return (-1,-1,-1,-1,-1)
    end
    #print("xoverall ")
    #println(xoverall)
    #print("yoverall ")
    #println(yoverall)

    #Calculate expected population size in next generation according to
    #Beverton-Holt
    Nexptd = 2*Rm*meanf*Nf*(1-(1-pmat)^Nm)/(1 + (Nf+Nm)/α);

    if Nexptd < 0 || isnan(Nexptd) == 1
        print("Nexptd = ")
        print(Nexptd)
        print("mean female fitness = ")
        print(meanf)
        print("Number of females = ")
        print(Nf)
        print("Number of males = ")
        print(Nm)
        print("α = ")
        print(α)
    end
    #Nnew = rand(Poisson(Nexptd));

    #Now sample according to Wright-Fisher
    if Nexptd > 10^15
        #Check with supervisor!!!
        Nmalenew = round(Int128,rand(Normal(ψ*Nexptd,sqrt(ψ*Nexptd))));
        Nfnew = round(Int128,rand(Normal((1-ψ)*Nexptd,sqrt((1-ψ)*Nexptd))));
    else
        Nmalenew = pois_rand(ψ*Nexptd);
        Nfnew = pois_rand((1-ψ)*Nexptd);
    end
    Nnew = Nmalenew + Nfnew;
    #Nmalenew = rand(Binomial(Nnew,ψ));
    #Nfnew = Nnew - Nmalenew;
    ###Alternatively, Nmalenew~Poi(ψ*Nexptd) and Nfnew~Poi((1-ψ)*Nexptd);
    if Gauss == 0
        Nfnew = round(Int, Nfnew);
        Nmalenew = round(Int, Nmalenew);
        xnumnew = rand(Multinomial(Nfnew,xoverall));
        ynumnew = rand(Multinomial(Nmalenew,yoverall));
    else
        xnumnew = GaussPoissonHybrid_mnrnd(Nfnew,xoverall);
        ynumnew = GaussPoissonHybrid_mnrnd(Nmalenew,yoverall);
    end
#nope!!
    if Nfnew > 0
        xnew = xnumnew/Nfnew;
    else
        xnew = zeros(na);
    end

    if Nmalenew > 0
        ynew = ynumnew/Nmalenew;
    else
        ynew = zeros(na);
    end

    return (xnew,ynew,Nfnew,Nmalenew,Nnew)
end
