include("functions.jl")

function test_ghost_points()
    """ Unit test for testing the ghost points function """

    matrix = [2 1 4 5;
              4 1 6 8;
              1 7 1 7;
              8 6 3 8]

    matrix_with_ghostpoints = [0 8 6 3 8 0;
                               5 2 1 4 5 2;
                               8 4 1 6 8 4;
                               7 1 7 1 7 1;
                               8 8 6 3 8 8
                               0 2 1 4 5 0]

    """ add ghost points with the written functions """
    matrix_ghostpoints_func = create_ghost_points(matrix)

    """ test equality """
    @test isequal(matrix_ghostpoints_func, matrix_with_ghostpoints)
end


function test_init()
    """
    Function for testing the init function
    against a 4x4 lattice with all spins up
     """

    """ Known values """
    L = 4
    E_known = -32
    M_known = 16

    spin_matrix_known = 1 .+ zeros(Int8, 4, 4)

    """ Compute spin matrix, energy and magnetization from init function """
    lattice, E, M = init(4, 1.0, true) # initialize as ordered

    """ extract inner matrix without ghost points """
    spin_matrix = zeros(Int8, 4, 4)
    for i = 2:L+1
        spin_matrix[i-1, :] = lattice[i,2:end-1]
    end

    @testset begin
        @test isequal(spin_matrix_known, spin_matrix)
        @test isequal(E_known, E)
        @test isequal(M_known, M)
    end

end

"""
julia> test_ghost_points()
Test Passed


julia> test_init()
Test Summary: | Pass  Total
test set      |    3      3
"""
