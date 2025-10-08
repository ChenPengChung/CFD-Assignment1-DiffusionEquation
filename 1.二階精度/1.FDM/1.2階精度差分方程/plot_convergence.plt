# Gnuplot script for grid convergence analysis
set terminal png enhanced size 800,600
set output 'grid_convergence.png'
set title 'Grid Convergence Analysis - L1 Error vs Grid Spacing'
set xlabel 'log(dx)'
set ylabel 'log(L1 Error)'
set grid
set key left top

# �u�ʦ^�k�u: y = 1.92982 * x + -2.54486
# �z��2����׽u�q�L�Ѧ��I (-2.99573, -8.29273)
f(x) = 1.92982 * x + -2.54486
g(x) = 2.0 * (x - -2.99573) + -8.29273

plot 'grid_convergence_data.dat' using 3:5 with linespoints pt 7 ps 1.5 lw 2 title sprintf('Computed (slope = %.2f)', 1.92982), \
     f(x) with lines lw 2 lc rgb 'red' title sprintf('Linear Fit (slope = %.2f)', 1.92982), \
     g(x) with lines lw 2 lc rgb 'green' dashtype 2 title '2nd Order Theory (slope = 2.0)'
