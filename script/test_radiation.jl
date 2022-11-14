include("../src/parameters.jl")
include("../src/setup.jl")
include("../src/drive.jl")
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



