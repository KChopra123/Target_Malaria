function unifmignewfrequency(Rm,α,xold,yold,μ,Wf,Wm,na,ψ,mutmatrix,κ,Nf,Nm,Gauss,pmig)
    # Updates haplotype frequency in next generation for each (square) deme
    # Reproduction occurs in each deme, and then the progeny migrate

    n = length(Nm);
    #println(1)

    xnew = zeros(n,na);
    ynew = zeros(n,na);
    Nfnew = zeros(Int128,n);
    Nmalenew = zeros(Int128,n);
    Nnew = zeros(Int128,n);
    #println(2)

    for iii = 1:n
            (xnew[iii,:],ynew[iii,:],Nfnew[iii],Nmalenew[iii],Nnew[iii]) = newfrequency(Rm,α,xold[iii,:],yold[iii,:],μ,Wf,Wm,na,ψ,mutmatrix,κ,Nf[iii],Nm[iii],Gauss);
    end
    #print(3)

    haplotype_number_x_new = zeros(Int128,(n,na));
    haplotype_number_y_new = zeros(Int128,(n,na));
    #print(4)

    for iii = 1:n
            for kkkk = 1:na
                # Check the following few lines!!!
                if isnan(xnew[iii,kkkk]*Nfnew[iii]) == 1
                    println(xold[iii,:])
                    println(yold[iii,:])
                    println(Nf[iii])
                    println(Nm[iii])
                    print(xnew[iii,:])
                    println(Nfnew[iii])
                    print(ynew[iii,:])
                    println(Nmalenew[iii])
                    #print("iii=")
                    #println(iii)
                    #print("jjj=")
                    #println(jjj)
                    #print("allele=")
                    #print(kkkk)
                    haplotype_number_x_new[iii,kkkk] = 0;
                    return (55,iii,kkkk)
                end
                if isnan(ynew[iii,kkkk]*Nmalenew[iii]) == 1
                    println(xold[iii,:])
                    println(yold[iii,:])
                    println(Nf[iii])
                    println(Nm[iii])
                    print(xnew[iii,:])
                    println(Nfnew[iii])
                    print(ynew[iii,:])
                    println(Nmalenew[iii])
                    println("allele=")
                    println(kkkk)
                    haplotype_number_y_new[iii,kkkk] = 0;
                    return (56,iii,kkkk)
                end
                haplotype_number_x_new[iii,kkkk] = round(Int128,xnew[iii,kkkk]*Nfnew[iii]);
                haplotype_number_y_new[iii,kkkk] = round(Int128,ynew[iii,kkkk]*Nmalenew[iii]);
            end
    end
    #print(5)

    # Now migrate the haplotypes independently
    B_x = zeros(Int128,(n,na));
    B_y = zeros(Int128,(n,na));
    #println(6)
    for kkkk = 1:na
        A_x = haplotype_number_x_new[:,kkkk];
        A_y = haplotype_number_y_new[:,kkkk];
        #println(6.5)

        B_x[:,kkkk] = unif_migration(A_x,pmig,Gauss);
        B_y[:,kkkk] = unif_migration(A_y,pmig,Gauss);
    end
    #print(7)

    # Finally, calculate new population sizes and haplotype frequencies
    Nf_post_mig = zeros(Int128,n);
    Nmale_post_mig = zeros(Int128,n);
    N_post_mig = zeros(Int128,n)
    x_post_mig = zeros(n,na);
    y_post_mig = zeros(n,na);
    #print(8)

    for iii = 1:n
            Nf_post_mig[iii] = sum(B_x[iii,:]);
            Nmale_post_mig[iii] = sum(B_y[iii,:]);
            N_post_mig[iii] = Nf_post_mig[iii] + Nmale_post_mig[iii];
            if Nf_post_mig[iii] > 0
                x_post_mig[iii,:] = B_x[iii,:]./Nf_post_mig[iii];
            else
                x_post_mig[iii,:] = B_x[iii,:];
            end

            if Nmale_post_mig[iii] > 0
                y_post_mig[iii,:] = B_y[iii,:]./Nmale_post_mig[iii];
            else
                y_post_mig[iii,:] = B_y[iii,:];
            end
    end
    #print(9)

    return (x_post_mig,y_post_mig,Nf_post_mig,Nmale_post_mig,N_post_mig,B_x,B_y)

    end
