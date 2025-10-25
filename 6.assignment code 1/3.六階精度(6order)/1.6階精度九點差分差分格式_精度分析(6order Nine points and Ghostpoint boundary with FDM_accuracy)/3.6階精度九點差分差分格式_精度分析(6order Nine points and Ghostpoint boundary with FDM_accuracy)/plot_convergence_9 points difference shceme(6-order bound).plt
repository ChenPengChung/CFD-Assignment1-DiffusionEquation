set terminal png enhanced size 800,600
set output 'grid_convergence_9 points difference shceme(6-order bound).png'
set title '9 points difference scheme (6-order bound)'
set xlabel 'log(dx)'
set ylabel 'log(L1 Error)'
set grid
set key left top
f(x) = 5.24533 * x + -6.10883
g(x) = 6.0 * (x - -2.99573) + -21.9187
plot 'grid_convergence_9 points difference shceme(6-order bound).dat' using 3:5 with linespoints pt 7 ps 1.5 lw 2 title sprintf('Improved (slope = %.2f)', 5.24533), \
     f(x) with lines lw 2 lc rgb 'red' title sprintf('Linear Fit (slope = %.2f)', 5.24533), \
     g(x) with lines lw 2 lc rgb 'green' dashtype 2 title '9 points difference shceme(6-order bound)(slope = 6.0)'
