# Gnuplot script for grid convergence analysis
set terminal png enhanced size 800,600
set output '1.2orderCD_FDM.png'
set title '1.2orderCD_FDM- L1 Error vs Grid Spacing'
set xlabel 'log(dx)'
set ylabel 'log(L1 Error)'
set grid
set key left top
f1(x) = 1.92982 * x + -2.54486
plot '1.2orderCD_FDM.dat' using 3:5  pt 0 ps 1.5 lc 'green' title sprintf('2orderCD_FDM (slope = %.2f)', 1.92982), \
     f1(x) with lines lw 2 lc 'green' title sprintf('2orderCD_FDM linear (slope = %.2f)', 1.92982), \
