#INDEX - Use ctrl+F to browse faster
#(1) - BINARY OPERATIONS
#(2) - SPIN OPERATIONS
##(2.1) - SPIN SITE OPERATIONS
##(2.2) - SYMMETRIES
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
                H[l,l]=H[l,l] - 1
            else
                H[l,l]=H[l,l] +1
            end             
        end
    end
    return H
end

##(2.2) SYMMETRIES

#S_z conservation - fixed spin up subspaces

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


function sz_subspace_spin_interactions(sites,n_partic,neig,BC)
    states, flag = sz_subspace(sites,n_partic)
    dim = convert(Int64,factorial(sites)/(factorial(sites-n_partic)*factorial(n_partic)))
    Sx = zeros(ComplexF64,dim,dim)
    Sy = zeros(ComplexF64,dim,dim)
    Sz = zeros(ComplexF64,dim,dim)
    Cx = 1
    Cy = 1
    Cz = 1
    t1 = 0
    kk = 0
    println(flag)
    println(typeof(flag))
    for i=1:dim
        for n=0:sites-2
            if ((n<=sites-1-neig) | (BC=="perodic"))
                stepi = ibits(states[i],n,1) + ibits(states[i], (n+neig)%(sites),1)
                println("stepi = $stepi")
                if (stepi == 1)
                    kk = ibits(states[i],n,1) + ibits(states[i],(n+neig)%(sites),1)
                    #t1 = (states[i])⊻(set_bit(0,n,true))
                    #println("first t1 = $t1")
                    t1 = ((states[i])⊻(set_bit(0,n,true)))⊻(set_bit(0,(n+neig)%(sites), true))
                    t1 = flag[t1+1]
                    #t1 = flag[t1]
                    #println(typeof(flag))
                    println("t1 = $t1")
                    Sy[i,t1]+=-Cy * (-1)^(kk)
                    Sx[i,t1]+= Cx                    
                end #if1
                Sz[i,i]+= Cz * (-1)^(ibits(states[i],n,1) + ibits(states[i],n+neig % sites,1))
            end #if2
        end #for2
    end #for
    return Sx + Sy + Sz
end #function

function cabeza_subspace(sites,n_partic)
    states, flag = sz_subspace(sites,n_partic)
    return S_zzij()
end

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

