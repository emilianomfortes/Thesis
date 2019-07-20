#INDEX - Use ctrl+F to browse faster
#(1) - BINARY OPERATIONS
#(2) - SPIN OPERATIONS
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

#----------------  (2) SPIN OPERATIONS  ----------------#

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

