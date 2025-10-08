set terminal png enhanced size 800,600
set output 'grid_convergence_9 points difference shceme(4-order bound).png'
set title 'Improved Grid Convergence Analysis: L1 Error vs Grid Spacing'
set xlabel 'log(dx)'
set ylabel 'log(L1 Error)'
set grid
set key left top
f(x) = 3.97604 * x + 0.234779
g(x) = 4.0 * (x - -2.99573) + -11.5937
plot 'grid_convergence_9 points difference shceme(4-order bound).dat' using 3:5 with linespoints pt 7 ps 1.5 lw 2 title sprintf('Improved (slope = %.2f)', 3.97604), \
     f(x) with lines lw 2 lc rgb 'red' title sprintf('Linear Fit (slope = %.2f)', 3.97604), \
     g(x) with lines lw 2 lc rgb 'green' dashtype 2 title '9 points difference shceme(4-order bound)(slope = 4.0)'
