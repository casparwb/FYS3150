The program is designed to the run in the Julia REPL terminal. To do this, the file containing the functions has to be included. 
This is done as follows in the REPL terminal:

julia> include("main.jl")

The functions can then be called directly in the terminal. For example:

julia> include("main.jl")
julia> CPU_time(100)

Produces the following output:

"""General algorithm time with 100 grid points: 0.002180 ms
Special algorithm time with 100 grid points: 0.000960 ms
Relative difference: 2.269651

(2.18e-6, 9.605e-7)"""