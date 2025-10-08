# Gnuplot script for grid convergence analysis
set terminal png enhanced size 800,600
set output '2.2orderCD_FVM.png'
set title '2.2orderCD_FVM - L1 Error vs Grid Spacing'
set xlabel 'log(dx)'
set ylabel 'log(L1 Error)'
set grid
set key left top
f2(x) = 1.9841 * x + -2.10977
plot '2.2orderCD_FVM.dat' using 3:5  pt 14 ps 1.5 lc 'orange' title sprintf('2orderCD_FVM (slope = %.2f)', 1.9841), \
     f2(x) with lines lw 2 lc 'orange' title sprintf('2orderCD_FVM linear (slope = %.2f)', 1.9841), \
     
