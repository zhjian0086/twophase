#bin/bash

#gnome-terminal -- bash -c "cd /mnt/HD01/Hossain/Code;
#cp HossainMPI /mnt/HD01/Hossain/R1000W2353/H001;
#cp HossainMPI /mnt/HD01/Hossain/R1000W2353/H002;
#cp HossainMPI /mnt/HD01/Hossain/R1000W2353/H005;
#cp HossainMPI /mnt/HD01/Hossain/R1000W2353/H010;
#cp HossainMPI /mnt/HD01/Hossain/R1000W2353/H020;
#cp HossainMPI /mnt/HD01/Hossain/R1000W2353/H050;
#cp HossainMPI /mnt/HD01/Hossain/R1000W2353/H100;
#cp HossainMPI /mnt/HD01/Hossain/R1000W2353/H200;"

gnome-terminal -- bash -c "cd /mnt/HOSSAIN/Pool-Depth/R1000W300/H2.00; mpiexec -n 5 ./DLnewmpi n6 x14 r1000.0 w300.0 h2.0;"
gnome-terminal -- bash -c "cd /mnt/HOSSAIN/Pool-Depth/R1000W300/H1.00; mpiexec -n 5 ./DLnewmpi n6 x14 r1000.0 w300.0 h1.0;"
gnome-terminal -- bash -c "cd /mnt/HOSSAIN/Pool-Depth/R1000W300/H0.50; mpiexec -n 5 ./DLnewmpi n6 x14 r1000.0 w300.0 h0.50;"
gnome-terminal -- bash -c "cd /mnt/HOSSAIN/Pool-Depth/R1000W300/H0.20; mpiexec -n 5 ./DLnewmpi n6 x14 r1000.0 w300.0 h0.20;"
gnome-terminal -- bash -c "cd /mnt/HOSSAIN/Pool-Depth/R1000W300/H0.10; mpiexec -n 5 ./DLnewmpi n6 x14 r1000.0 w300.0 h0.10;"

#gnome-terminal -x bash -c "cd /mnt/HD01/HOSSAIN/PoolDepthEffect/TestCases/R2000FS001/H001; ./DLpost 6 13 0.0 0.001 0.95; /usr/local/tecplot360ex/bin/tec360 -b -p CC-Interface.mcr; ./JetRecognition 6 13 0.0 0.001 0.95"
#gnome-terminal -x bash -c "cd /mnt/HD01/HOSSAIN/PoolDepthEffect/TestCases/R2000FS001/H002; ./DLpost 6 14 0.0 0.001 0.95; /usr/local/tecplot360ex/bin/tec360 -b -p CC-Interface.mcr; ./JetRecognition 6 14 0.0 0.001 0.95"
#gnome-terminal -x bash -c "cd /mnt/HD01/HOSSAIN/PoolDepthEffect/TestCases/R2000FS001/H005; ./DLpost 6 13 0.0 0.001 0.95; /usr/local/tecplot360ex/bin/tec360 -b -p CC-Interface.mcr; ./JetRecognition 6 13 0.0 0.001 0.95"
#gnome-terminal -x bash -c "cd /mnt/HD01/HOSSAIN/PoolDepthEffect/TestCases/R2000FS001/H010; ./DLpost 6 13 0.0 0.001 0.95; /usr/local/tecplot360ex/bin/tec360 -b -p CC-Interface.mcr; ./JetRecognition 6 13 0.0 0.001 0.95"
#gnome-terminal -x bash -c "cd /mnt/HD01/HOSSAIN/PoolDepthEffect/TestCases/R2000FS001/H020; ./DLpost 6 13 0.0 0.001 0.95; /usr/local/tecplot360ex/bin/tec360 -b -p CC-Interface.mcr; ./JetRecognition 6 13 0.0 0.001 0.95"
#gnome-terminal -x bash -c "cd /mnt/HD01/HOSSAIN/PoolDepthEffect/TestCases/R2000FS001/H050; ./DLpost 6 13 0.0 0.001 0.95; /usr/local/tecplot360ex/bin/tec360 -b -p CC-Interface.mcr; ./JetRecognition 6 13 0.0 0.001 0.95"
#gnome-terminal -x bash -c "cd /mnt/HD01/HOSSAIN/PoolDepthEffect/TestCases/R2000FS001/H100; ./DLpost 6 13 0.0 0.001 0.95; /usr/local/tecplot360ex/bin/tec360 -b -p CC-Interface.mcr; ./JetRecognition 6 13 0.0 0.001 0.95"
#gnome-terminal -x bash -c "cd /mnt/HD01/HOSSAIN/PoolDepthEffect/TestCases/R2000FS001/H200; ./DLpost 6 13 0.0 0.001 0.95; /usr/local/tecplot360ex/bin/tec360 -b -p CC-Interface.mcr; ./JetRecognition 6 13 0.0 0.001 0.95"

#gnome-terminal -- bash -c "cd /media/hchizari/HD01/Hossain/SaveFile001-R1000W300H001; ./JetRecognition 6 13 0.0 0.001 0.95"
#gnome-terminal -- bash -c "cd /media/hchizari/HD01/Hossain/SaveFile001-R1000W300H002; ./JetRecognition 6 13 0.0 0.001 0.95"
#gnome-terminal -- bash -c "cd /media/hchizari/HD01/Hossain/SaveFile001-R1000W300H005; ./JetRecognition 6 13 0.0 0.001 0.95"
#gnome-terminal -- bash -c "cd /media/hchizari/HD01/Hossain/SaveFile001-R1000W300H010; ./JetRecognition 6 13 0.0 0.001 0.95"
#gnome-terminal -- bash -c "cd /media/hchizari/HD01/Hossain/SaveFile001-R1000W300H020; ./JetRecognition 6 13 0.0 0.001 0.95"
#gnome-terminal -- bash -c "cd /media/hchizari/HD01/Hossain/SaveFile001-R1000W300H050; ./JetRecognition 6 14 0.0 0.001 0.95"
#gnome-terminal -- bash -c "cd /media/hchizari/HD01/Hossain/SaveFile001-R1000W300H100; ./JetRecognition 6 13 0.0 0.001 0.95"
#gnome-terminal -- bash -c "cd /media/hchizari/HD01/Hossain/SaveFile001-R1000W300H200; ./JetRecognition 6 13 0.0 0.001 0.95"
