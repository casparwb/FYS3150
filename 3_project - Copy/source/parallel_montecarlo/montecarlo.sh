for i in 1000 10000 100000 1000000 10000000 100000000; do
./mc_serial.exe $i serial_n=$i.txt
mpirun -np 4 ./mc_MPI.exe $i parallel_n=$i.txt
done
