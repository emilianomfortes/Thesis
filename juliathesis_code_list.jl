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
