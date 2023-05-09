function sqmignewfrequency_allee(Rm,α,xold,yold,μ,Wf,Wm,na,ψ,mutmatrix,κ,Nf,Nm,Gauss,pmig,pmat)
    # Updates haplotype frequency in next generation for each (square) deme
    # Reproduction occurs in each deme, and then the progeny migrate

    (acrm,upn) = size(Nm);
    #println(1)

    xnew = zeros(acrm,upn,na);
    ynew = zeros(acrm,upn,na);
    Nfnew = zeros(Int128,(acrm,upn));
    Nmalenew = zeros(Int128,(acrm,upn));
    Nnew = zeros(Int128,(acrm,upn));
    #println(2)

    for iii = 1:acrm
        for jjj = 1:upn
            (xnew[iii,jjj,:],ynew[iii,jjj,:],Nfnew[iii,jjj],Nmalenew[iii,jjj],Nnew[iii,jjj]) = newfrequency_allee(Rm,α,xold[iii,jjj,:],yold[iii,jjj,:],μ,Wf,Wm,na,ψ,mutmatrix,κ,Nf[iii,jjj],Nm[iii,jjj],pmat,Gauss);
        end
    end
    #print(3)

    haplotype_number_x_new = zeros(Int128,(acrm,upn,na));
    haplotype_number_y_new = zeros(Int128,(acrm,upn,na));
    #print(4)

    for iii = 1:acrm
        for jjj = 1:upn
            for kkkk = 1:na
                # Check the following few lines!!!
                if isnan(xnew[iii,jjj,kkkk]*Nfnew[iii,jjj]) == 1
                    println(xold[iii,jjj,:])
                    println(yold[iii,jjj,:])
                    println(Nf[iii,jjj])
                    println(Nm[iii,jjj])
                    print(xnew[iii,jjj,:])
                    println(Nfnew[iii,jjj])
                    print(ynew[iii,jjj,:])
                    println(Nmalenew[iii,jjj])
                    #print("iii=")
                    #println(iii)
                    #print("jjj=")
                    #println(jjj)
                    #print("allele=")
                    #print(kkkk)
                    haplotype_number_x_new[iii,jjj,kkkk] = 0;
                    return (55,iii,jjj,kkkk)
                end
                if isnan(ynew[iii,jjj,kkkk]*Nmalenew[iii,jjj]) == 1
                    println(xold[iii,jjj,:])
                    println(yold[iii,jjj,:])
                    println(Nf[iii,jjj])
                    println(Nm[iii,jjj])
                    print(xnew[iii,jjj,:])
                    println(Nfnew[iii,jjj])
                    print(ynew[iii,jjj,:])
                    println(Nmalenew[iii,jjj])
                    println("allele=")
                    println(kkkk)
                    haplotype_number_y_new[iii,jjj,kkkk] = 0;
                    return (56,iii,jjj,kkkk)
                end
                haplotype_number_x_new[iii,jjj,kkkk] = round(Int128,xnew[iii,jjj,kkkk]*Nfnew[iii,jjj]);
                haplotype_number_y_new[iii,jjj,kkkk] = round(Int128,ynew[iii,jjj,kkkk]*Nmalenew[iii,jjj]);
            end
        end
    end
    #print(5)

    # Now migrate the haplotypes independently
    B_x = zeros(Int128,(acrm,upn,na));
    B_y = zeros(Int128,(acrm,upn,na));
    #println(6)
    for kkkk = 1:na
        A_x = haplotype_number_x_new[:,:,kkkk];
        A_y = haplotype_number_y_new[:,:,kkkk];
        #println(6.5)

        B_x[:,:,kkkk] = square_deme_migration(A_x,pmig,Gauss);
        B_y[:,:,kkkk] = square_deme_migration(A_y,pmig,Gauss);
    end
    #print(7)

    # Finally, calculate new population sizes and haplotype frequencies
    Nf_post_mig = zeros(Int128,(acrm,upn));
    Nmale_post_mig = zeros(Int128,(acrm,upn));
    N_post_mig = zeros(Int128,(acrm,upn))
    x_post_mig = zeros(acrm,upn,na);
    y_post_mig = zeros(acrm,upn,na);
    #print(8)

    for iii = 1:acrm
        for jjj = 1:upn
            Nf_post_mig[iii,jjj] = sum(B_x[iii,jjj,:]);
            Nmale_post_mig[iii,jjj] = sum(B_y[iii,jjj,:]);
            N_post_mig[iii,jjj] = Nf_post_mig[iii,jjj] + Nmale_post_mig[iii,jjj];
            if Nf_post_mig[iii,jjj] > 0
                x_post_mig[iii,jjj,:] = B_x[iii,jjj,:]./Nf_post_mig[iii,jjj];
            else
                x_post_mig[iii,jjj,:] = B_x[iii,jjj,:];
            end

            if Nmale_post_mig[iii,jjj] > 0
                y_post_mig[iii,jjj,:] = B_y[iii,jjj,:]./Nmale_post_mig[iii,jjj];
            else
                y_post_mig[iii,jjj,:] = B_y[iii,jjj,:];
            end
        end
    end
    #print(9)

    return (x_post_mig,y_post_mig,Nf_post_mig,Nmale_post_mig,N_post_mig,B_x,B_y)

    end
