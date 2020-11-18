two = BigInt(2)
ten = BigInt(10)

for i in 1:9
    a = BigInt(i)

    for j in 2:50
        a += i*two^(j-1)*ten^(j-1)

        at = a % ten^j
        s = string(at)
        si = s[end]*s[1:end-1]
        s2 = string(2at)

        if si == s2
            println(s)
        end
    end
end
