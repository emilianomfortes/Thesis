function btest(i,pos)
    a = i & (1 << pos)
    if a >= 1
        a = 1
    end
    return Core.Bool(a)
end

function set_bit(v,index,x)
    mask = 1 << index
    v &= ~mask
    if x==true
        v |= mask
    end
    return v
end
