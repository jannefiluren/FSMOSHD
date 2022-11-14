tmp = readline("../fortran/data/input_fake_5wj.txt")
tmp = parse.(Float64,split(tmp))

year = tmp[1]
month = tmp[2]
day = tmp[3]
hour = tmp[4]
Sdir[:,:] .= tmp[5]
Sdif[:,:] .= tmp[6]
LW[:,:] .= tmp[7]
Sf[:,:] .= tmp[8]
Rf[:,:] .= tmp[9]
Ta[:,:] .= tmp[10]
RH[:,:] .= tmp[11]
Ua[:,:] .= tmp[12]
Ps[:,:] .= tmp[13]
Sf24h[:,:] .= tmp[14]
