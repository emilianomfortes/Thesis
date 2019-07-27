using LinearAlgebra
using StatsBase
using Distribution
using LsqFit
#INDEX - Use ctrl+F to browse faster
#(1) - BINARY OPERATIONS
#(2) - SPIN OPERATIONS
##(2.1) - SPIN SITE OPERATIONS
##(2.2) - CHAIN MODELS
##(2.3) - SYMMETRIES
###(2.3.1) - S_z CONSERVATION
#(3) - OUT-OF-TIME-ORDERED CORRELATORS
##(3.1) - WITH TEMPERATURE CHOICES (IF T=INFTY USE INFTY TEMPERATURE OTOCS FOR EXTRA SPEED)
##(3.2) - INFTY TEMPERATURE OTOCS
#(4) - CHAOS OPERATIONS
#(5) - MISCELLANEA CODE

#----------------  (1) BINARY OPERATIONS  ----------------#
#Translations from fortran/python --> julia

#btest (https://gnu.huihoo.org/gcc/gcc-7.1.0/gfortran/BTEST.html)
function btest(i,pos)
    a = i & (1 << pos)
    if a >= 1
        a = 1
    end
    return Core.Bool(a)
end

#ibclr/ibset (https://gcc.gnu.org/onlinedocs/gcc-4.6.1/gfortran/IBCLR.html / http://www.lahey.com/docs/lfpro78help/F95ARIBSETFn.htm)
function set_bit(v,index,x)
    mask = 1 << index
    v &= ~mask
    if x==true
        v |= mask
    end
    return v
end


struct BitStrFun2{L} # L is the parameter which gets the string length
    value::Int           
    function BitStrFun2(str::String) where {L}
        value = parse(Int, string("0b", str))
        nbits = length(str)
        return new{nbits}(value)
     end
end

#ibits fortran function
function ibits(num,pos,lens)
    binary = Base.reverse(Base.bitstring(num))
    return BitStrFun2(Base.reverse(binary[pos+1:pos+1+lens-1])).value
end

#ieor - I LEAVE IT HERE BECAUSE IT TOOK ME AN HOUR TO CODE IT, BUT ITS AN IMPLICIT FUNCTION " ⊻ " (IE 1⊻4 = 5)
function ieor(num1,num2)
    a = length(Base.digits(num1, base=2))
    b = length(Base.digits(num2, base=2))
    lenti = max(a,b)
    print(" len(num1) = $a and len(num2) = $b")
    if lenti == a
        zz = zeros(Int,lenti)
        zz[1:b] = Base.digits(num2,base=2)
        println(zz)
        numeronuevo = abs.(zz-Base.digits(num1,base=2))
    else
        zz = zeros(Int,lenti)
        zz[1:a] = Base.digits(num1,base=2)
        println(zz)
        numeronuevo = abs.(zz-Base.digits(num2,base=2))
    end
    println(numeronuevo)
    println(join(numeronuevo))
    Base.parse(Int,join(Base.reverse(numeronuevo)),base=2)
end


#----------------  (2) SPIN OPERATIONS  ----------------#

##(2.1) SPIN SITE OPERATIONS

# Pauli at site operators
# Pauli at site operator (x,i) - The opt version is faster than the non opt version if sites>10
function S_xi_opt(pos_i,sites)
    dim = 2^sites
    S = zeros(ComplexF64,dim,dim)
    for i=0:dim-1
        t1 = (i)⊻(set_bit(0,pos_i,true))
        S[i+1,t1+1]+=1
    end
    return S
end #function

function S_xi(pos_i,sites)
    dim = 2^sites
    S = zeros(ComplexF64,dim,dim)
    estados2 = zeros(Int64,dim)
    for i=0:dim-1
        if btest(i,pos_i) == true
            estados2[i+1] = set_bit(i,pos_i,0)
        else
            estados2[i+1] = set_bit(i,pos_i,true)
        end
    end
    for i=0:dim-1
        for j=0:dim-1
            if i == estados2[j+1]
                S[i+1,j+1] = S[i+1,j+1]+1
            end
        end
    end
    return S
end

# Pauli at site operators (y,i) - The opt version is faster than the non opt version if sites>10
function S_yi_opt(pos_i,sites)
    dim = 2^sites
    S = zeros(ComplexF64,dim,dim)
    for i=0:dim-1
        t1 = (i)⊻(set_bit(0,pos_i,true))
        S[i+1,t1+1]+= -1im*((-1)^(ibits(i,pos_i,1)))
    end
    return S
end #function

function S_yi(pos_i,sites)
    dim = 2^sites
    S = zeros(ComplexF64,dim,dim)
    estados2 = zeros(Int64,dim)
    a = zeros(ComplexF64,dim)
    for i=0:dim-1
        if btest(i,pos_i) == true
            estados2[i+1] = set_bit(i,pos_i,0)
            a[i+1] = 1im
        else
            estados2[i+1] = set_bit(i,pos_i,true)
            a[i+1] = -1im
        end
    end
    for i=0:dim-1
        for j=0:dim-1
            if i == estados2[j+1]
                S[i+1,j+1] = S[i+1,j+1]+a[i+1]
            end
        end
    end
    return S
end

# Pauli at site operators (z,i) - The opt version is faster than the non opt version if sites>10
function S_zi_opt(pos_i,sites)
    dim = 2^sites
    S = zeros(ComplexF64,dim,dim)
    for i=0:dim-1
        S[i+1,i+1]+= (-1)^(ibits(i,pos_i,1))
    end
    return S
end #function

function S_zi(pos_i,sites)
    dim = 2^sites
    S = zeros(ComplexF64,dim,dim)
    for i=0:dim-1
        if btest(i,pos_i) == false
            S[i+1,i+1] = 1
        else
            S[i+1,i+1] = -1
        end
    end
    return S
end

#Entire spin direction operators

function S_x(sites)
    dim = 2^sites
    H = zeros(ComplexF64,dim,dim)
    for l=1:dim
        for j=0:sites-1
            if btest(l-1,j) == true
                hh = set_bit(l-1,j,0) + 1
            else
                hh = set_bit(l-1,j,true) + 1
            end             
            H[l,hh] = H[l,hh] + 1
        end
    end
    return H
end

function S_y(sites)
    dim = 2^sites
    H = zeros(ComplexF64,dim,dim)
    for l=1:dim
        for j=0:sites-1
            if btest(l-1,j) == true
                hh = set_bit(l-1,j,0) + 1
                H[l,hh] = H[l,hh] + 1im
            else
                hh = set_bit(l-1,j,true) + 1
                H[l,hh] = H[l,hh] - 1im
            end             
        end
    end
    return H
end

function S_z(sites)
    dim = 2^sites
    H = zeros(ComplexF64,dim,dim)
    for l=1:dim
        for j=0:sites-1
            if btest(l-1,j) == true
                H[l,l]+=- 1
            else
                H[l,l]+= +1
            end             
        end
    end
    return H
end

# Neighbour interactionss xx,yy,zz all in one for speed
function spin_interactions(sites,neig,BC,Cx,Cy,Cz)
    dim = 2^sites
    Sx = zeros(ComplexF64,dim,dim)
    Sy = zeros(ComplexF64,dim,dim)
    Sz = zeros(ComplexF64,dim,dim)
    t1 = 0
    kk = 0
    for i=0:dim-1
        for n=0:sites-2
            if ((n<=sites-1-neig) | (BC=="perodic"))
                kk = ibits(i,n,1) + ibits(i,(n+neig)%(sites),1)
                t1 = ((i)⊻(set_bit(0,n,true)))⊻(set_bit(0,(n+neig)%(sites), true))
                println("t1 = $t1")
                Sy[i+1,t1+1]+=-Cy * (-1)^(kk)
                Sx[i+1,t1+1]+= Cx                    
                Sz[i+1,i+1]+= Cz * (-1)^(kk)
            end #if
        end #for2
    end #for
    return Sx + Sy + Sz
end #function-

## (2.2) CHAIN MODELS

# Ising model with tilted magnetic field
function Tilted_model(sites,BC,J,B,tita)
    H = .25*J*spin_interactions(sites,1,BC,0,0,1)
    H+= .5*B*(sin(tita)*S_x(sites) + cos(tita)*S_z(sites))
    Hh = Hermitian(H)
    e = eigvals(H)
    ev = eigvecs(H)
    return e, ev
end #function

# Heisenberg chain with random magnetic field interactions
function Heisenberg_unif_random(sites,BC,h)
    H = .25*spin_interactions(sites,1,BC,1,1,1)
    hh = rand(sites) * 2h - h
    for i=0:sites-1
        H += .5* hh[i+1] * S_zi(i,sites)
    end #for
    #H += .5*B*(sin(tita)*S_x(sites) + cos(tita)*S_z(sites))
    Hh = Hermitian(H)
    e = eigvals(H)
    ev = eigvecs(H)
    return e, ev
end #function

# Perturbed XXZ model
function Perturbed_XXZ(sites,BC,lambdi,a)
    H = .25*spin_interactions(sites,1,BC,1,1,a)
    H += .25*lambdi*spin_interactions(sites,2,BC,1,1,a)
    Hh = Hermitian(H)
    e = eigvals(H)
    ev = eigvecs(H)
    return e, ev
end #function


##(2.3) SYMMETRIES

###(2.3.1) S_z CONSERVATION - fixed spin up subspace

#S_z conservation system states
function sz_subspace(sites,n_partic)
    label_state = 0
    dim = convert(Int64,factorial(sites)/(factorial(sites-n_partic)*factorial(n_partic)))
    dim2 = 2^sites
    states = zeros(Int64,dim)
    flag = zeros(Int64,dim2)
    for i=0:dim2-1
        k=0
        j=0
        tusi=0
        while ((k<=n_partic) & (j<=sites-1))
            k+=ibits(i,j,1)
            j+=1
            #println("k, n_partic = $k $n_partic")
            #println("j, sites-1 = $j $(sites-1)")
            #println("$k")
        end #while1
        #println("--------")
        if k==n_partic
            label_state+=1
            states[label_state] = i
            flag[i+1] = label_state
        end #if
    end #for
    return states, flag
end #function

# Pauli at site operator (x,i)
function sz_subspace_S_xi_opt(pos_i,sites,n_partic)
    dim = convert(Int64,factorial(sites)/(factorial(sites-n_partic)*factorial(n_partic)))
    S = zeros(ComplexF64,dim,dim)
    states, flag = sz_subspace(sites,n_partic)
    for i=1:dim
            t1 = (states[i])⊻(set_bit(0,pos_i,true))
            S[i,i]+=1
    end #for  
    return S
end

function sz_subspace_S_xi(pos_i,sites,n_partic)
    dim = convert(Int64,factorial(sites)/(factorial(sites-n_partic)*factorial(n_partic)))
    S = zeros(ComplexF64,dim,dim)
    states, flag = sz_subspace(sites,n_partic)
    estados2 = zeros(Int64,dim)
    for i=0:dim-1
        if btest(states[i+1],pos_i) == true
            estados2[i+1] = set_bit(states[i+1],pos_i,0)
        else
            estados2[i+1] = set_bit(states[i+1],pos_i,true)
        end #if  
    end
    for i=1:dim
        for j=1:dim
            if states[i] == estados2[j]
                S[i,j]+=1
            end #if
        end#for2
    end #for
    return S
end

# Pauli at site operator (y,i)
function sz_subspace_S_yi_opt(pos_i,sites,n_partic)
    dim = convert(Int64,factorial(sites)/(factorial(sites-n_partic)*factorial(n_partic)))
    S = zeros(ComplexF64,dim,dim)
    states, flag = sz_subspace(sites,n_partic)
    for i=1:dim
            t1 = (states[i])⊻(set_bit(0,pos_i,true))
            S[i,t1+1]+= -1im*((-1)^(ibits(states[i],pos_i,1)))
    end #for  
    return S
end

function sz_subspace_S_yi(pos_i,sites,n_partic)
    dim = convert(Int64,factorial(sites)/(factorial(sites-n_partic)*factorial(n_partic)))
    S = zeros(ComplexF64,dim,dim)
    states, flag = sz_subspace(sites,n_partic)
    estados2 = zeros(Int64,dim)
    a = zeros(ComplexF64,dim)
    for i=1:dim
        if btest(states[i],pos_i) == true
            estados2[i] = set_bit(states[i],pos_i,0)
            a[i] = 1im
        else
            estados2[i] = set_bit(states[i],pos_i,true)
            a[i] = -1im
        end #if  
    end
    for i=1:dim
        for j=1:dim
            if states[i] == estados2[j]
                S[i,j]+=a[i]
            end #if
        end#for2
    end #for
    return S
end

# Pauli at site operator (z,i)
function sz_subspace_S_zi_opt(pos_i,sites,n_partic)
    dim = convert(Int64,factorial(sites)/(factorial(sites-n_partic)*factorial(n_partic)))
    S = zeros(ComplexF64,dim,dim)
    states, flag = sz_subspace(sites,n_partic)
    for i=1:dim
            S[i,i]+= (-1)^(ibits(states[i],pos_i,1))
    end #for  
    return S
end

function sz_subspace_S_zi(pos_i,sites,n_partic)
    dim = convert(Int64,factorial(sites)/(factorial(sites-n_partic)*factorial(n_partic)))
    S = zeros(ComplexF64,dim,dim)
    states, flag = sz_subspace(sites,n_partic)
    for i=0:dim-1
        if btest(states[i+1],pos_i) == false
            S[i+1,i+1] = 1
        else
            S[i+1,i+1] = -1
        end #if  
    end
    return S
end


# Neighbor (x,i) (x,j) interaction
function sz_subspace_S_xxij(pos_i,pos_j,sites,n_partic)
    dim = convert(Int64,factorial(sites)/(factorial(sites-n_partic)*factorial(n_partic)))
    S = zeros(ComplexF64,dim,dim)
    states, flag = sz_subspace(sites,n_partic)
    estados2 = zeros(Int64,dim)
    for i=0:dim-1
        if btest(states[i+1],pos_i) == true
            estados2[i+1] = set_bit(states[i+1],pos_i,0)
        else
            estados2[i+1] = set_bit(states[i+1],pos_i,true)
        end
    end
    for i=0:dim-1
        if btest(estados2[i+1],pos_j) == true
            estados2[i+1] = set_bit(estados2[i+1],pos_j,0)
        else
            estados2[i+1] = set_bit(estados2[i+1],pos_j,true)
        end
    end
    for i=0:dim-1
        for j=0:dim-1
            if states[i+1] == estados2[j+1]
                S[i+1,j+1] = S[i+1,j+1]+1
            end
        end
    end
    return S
end

# Neighbor (y,i) (y,j) interaction

function sz_subspace_S_yyij(pos_i,pos_j,sites,n_partic)
    dim = convert(Int64,factorial(sites)/(factorial(sites-n_partic)*factorial(n_partic)))
    S = zeros(ComplexF64,dim,dim)
    states, flag = sz_subspace(sites,n_partic)
    estados2 = zeros(Int64,dim)
    a = zeros(ComplexF64,dim)
    for i=0:dim-1
        if btest(states[i+1],pos_i) == true
            estados2[i+1] = set_bit(states[i+1],pos_i,0)
            a[i+1] = 1im
        else
            estados2[i+1] = set_bit(states[i+1],pos_i,true)
        a[i+1] = -1im
        end
    end
    for i=0:dim-1
        if btest(estados2[i+1],pos_j) == true
            estados2[i+1] = set_bit(estados2[i+1],pos_j,0)
            a[i+1] = a[i+1] * 1im
        else
            estados2[i+1] = set_bit(estados2[i+1],pos_j,true)
            a[i+1] = -a[i+1]*1im
        end
    end
    for i=0:dim-1
        for j=0:dim-1
            if states[i+1] == estados2[j+1]
                S[i+1,j+1] = S[i+1,j+1]+a[i+1]
            end
        end
    end
    return S
end

# Neighbor (z,i) (z,j) interaction

function sz_subspace_S_zzij(pos_i,pos_j,sites,n_partic)
    dim = convert(Int64,factorial(sites)/(factorial(sites-n_partic)*factorial(n_partic)))
    S = zeros(ComplexF64,dim,dim)
    states, flag = sz_subspace(sites,n_partic)
    for i=0:dim-1
        if btest(states[i+1],pos_i) == true
            if btest(states[i+1],pos_j) == true
                S[i+1,i+1]+= 1
            else
                S[i+1,i+1]+= -1
            end
        else
            if btest(states[i+1],pos_j) == true
                S[i+1,i+1]+= -1
            else
                S[i+1,i+1]+= 1
            end
        end
    end
    return S
end

# Neighbor Cx (x,i)(x,j) + Cy(y,i)(y,j) + Cx(z,i)(z,j) interaction

function sz_subspace_spin_interactions(sites,n_partic,neig,BC,Cx,Cy,Cz)
    states, flag = sz_subspace(sites,n_partic)
    dim = convert(Int64,factorial(sites)/(factorial(sites-n_partic)*factorial(n_partic)))
    Sx = zeros(ComplexF64,dim,dim)
    Sy = zeros(ComplexF64,dim,dim)
    Sz = zeros(ComplexF64,dim,dim)
    t1 = 0
    kk = 0
    for i=1:dim
        for n=0:sites-2
            if ((n<=sites-1-neig) | (BC=="perodic"))
                stepi = ibits(states[i],n,1) + ibits(states[i], (n+neig)%(sites),1)
                if (stepi == 1)
                    kk = ibits(states[i],n,1) + ibits(states[i],(n+neig)%(sites),1)
                    t1 = ((states[i])⊻(set_bit(0,n,true)))⊻(set_bit(0,(n+neig)%(sites), true))
                    t1 = flag[t1+1]
                    Sy[i,t1]+=-Cy * (-1)^(kk)
                    Sx[i,t1]+= Cx                    
                end #if1
                Sz[i,i]+= Cz * (-1)^(ibits(states[i],n,1) + ibits(states[i],n+neig % sites,1))
            end #if2
        end #for2
    end #for
    return Sx + Sy + Sz
end #function

# Heisenberg chain with random uniform magnetic field in the Sz_subspace
function sz_subspace_Heisenberg_unif_random(sites,n_partic,h)
    H = .25*sz_subspace_spin_interactions(sites,n_partic,1,"OPEN",1,1,1)
    hh = rand(sites) * 2h - h
    for i=0:sites-1
        H += .5* hh[i+1] * sz_subspace_S_zi(i,sites)
    end #for
    Hh = Hermitian(H)
    e = eigvals(H)
    ev = eigvecs(H)
    return e, ev
end #function

#Perturbed XXZ in the S_z subspace
function sz_subspace_Perturbed_XXZ(sites,n_partic,BC,lambdi,a)
    H = .25*sz_subspace_spin_interactions(sites,n_partic,1,BC,1,1,a)
    H += .25*lambdi*sz_subspace_spin_interactions(sites,n_partic,2,BC,1,1,a)
    Hh = Hermitian(H)
    e = eigvals(H)
    ev = eigvecs(H)
    return e, ev
end #function

#----------------  (4) OUT-OF-TIME-ORDERED CORRELATORS  ----------------#

function OTOCF_infty(V,W,ener,basis,N,dt,t0,ortho)
    dim = length(ener)
    if ortho == true
        basist = transpose(basis)
    else
        basist = inv(basis)
    end#if
    ope0 = basist*V
    ope0 = ope0*basis
    ope = basist*W
    ope = ope*basis
    otoc = zeros(Float32,N)
    mm = zeros(ComplexF64,dim,dim)
    mm1 = zeros(ComplexF64,dim,dim)
    U = zeros(ComplexF64,dim,dim)
    Udag = zeros(ComplexF64,dim,dim)
    tiempo = LinRange(t0,N*dt+t0,N)
    for ti=1:N
        for c1=1:dim
            U[c1,c1] = exp(-1im*tiempo[ti]*ener[c1])
            Udag[c1,c1] = exp(1im*tiempo[ti]*ener[c1])
        end#for2
        #println(U)
        mm = ope0*U
        mm1 = Udag*mm
        mm = mm1*ope
        mm = mm*mm
        #println(tr(mm))
        otoc[ti] = 1 - Base.real((tr(mm)))/dim
    end#for
    return otoc,tiempo
end #function

