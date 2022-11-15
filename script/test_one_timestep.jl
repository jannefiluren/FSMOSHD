include("../src/parameters.jl")
include("../src/setup.jl")
include("../src/drive.jl")
include("../src/qsat.jl")

### TEST RADIATION

println("TEST RADIATION")

include("../src/radiation.jl")

fortran = readlines("../fortran/data/test_radiation.txt")

alb_fortran = parse.(Float64,split(fortran[1]))
asrf_out_fortran = parse.(Float64,split(fortran[2]))
Sdirt_fortran = parse.(Float64,split(fortran[3]))
Sdift_fortran = parse.(Float64,split(fortran[4]))
SWveg_fortran = parse.(Float64,split(fortran[5]))
SWsrf_fortran = parse.(Float64,split(fortran[6]))
SWsci_fortran = parse.(Float64,split(fortran[7]))
LWt_fortran = parse.(Float64,split(fortran[8]))

println(maximum(abs.(alb_fortran - alb)))
println(maximum(abs.(asrf_out_fortran - asrf_out)))
println(maximum(abs.(Sdirt_fortran - Sdirt)))
println(maximum(abs.(Sdift_fortran - Sdift)))
println(maximum(abs.(SWveg_fortran - SWveg)))
println(maximum(abs.(SWsrf_fortran - SWsrf)))
println(maximum(abs.(SWsci_fortran - SWsci)))
println(maximum(abs.(LWt_fortran - LWt)))

### TEST THERMAL

println("TEST THERMAL")

include("../src/thermal.jl")

fortran = readlines("../fortran/data/test_thermal.txt")

Ds1_fortran = parse.(Float64,split(fortran[1]))
gs1_fortran = parse.(Float64,split(fortran[2]))
ks1_fortran = parse.(Float64,split(fortran[3]))
Ts1_fortran = parse.(Float64,split(fortran[4]))
Tveg0_fortran = parse.(Float64,split(fortran[5]))
csoil_fortran = parse.(Float64,split(fortran[6]))
ksnow_fortran = parse.(Float64,split(fortran[7]))
ksoil_fortran = parse.(Float64,split(fortran[8]))

println(maximum(abs.(Ds1_fortran - Ds1)))
println(maximum(abs.(gs1_fortran - gs1)))
println(maximum(abs.(ks1_fortran - ks1)))
println(maximum(abs.(Ts1_fortran - Ts1)))
println(maximum(abs.(Tveg0_fortran - Tveg0)))
println(maximum(abs.(csoil_fortran - csoil[:,1,1])))
println(maximum(abs.(ksnow_fortran - ksnow[:,1,1])))
println(maximum(abs.(ksoil_fortran - ksoil[:,1,1])))

### TEST SFEXCH

println("TEST SFEXCH")

include("../src/sfexch.jl")

fortran = readlines("../fortran/data/test_sfexch.txt")

KH_fortran = parse.(Float64,split(fortran[1]))
KHa_fortran = parse.(Float64,split(fortran[2]))
KHg_fortran = parse.(Float64,split(fortran[3]))
KHv_fortran = parse.(Float64,split(fortran[4]))
KWg_fortran = parse.(Float64,split(fortran[5]))
KWv_fortran = parse.(Float64,split(fortran[6]))
Usc_fortran = parse.(Float64,split(fortran[7]))

println(maximum(abs.(KH_fortran - KH)))
println(maximum(abs.(KHa_fortran - KHa)))
println(maximum(abs.(KHg_fortran - KHg)))
println(maximum(abs.(KHv_fortran - KHv)))
println(maximum(abs.(KWg_fortran - KWg)))
println(maximum(abs.(KWv_fortran - KWv)))
println(maximum(abs.(Usc_fortran - Usc)))
