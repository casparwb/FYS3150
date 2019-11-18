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

    matrix_ghostpoints_func = create_ghost_points(matrix)

    @assert isequal(matrix_ghostpoints_func, matrix_with_ghostpoints) " Ghost point function erorr"
end
