# Gnuplot script for grid convergence analysis
set terminal png enhanced size 800,600
set output 'grid_convergence.png'
set title 'Grid Convergence Analysis - L1 Error vs Grid Spacing'
set xlabel 'log(dx)'
set ylabel 'log(L1 Error)'
set grid
set key left top

# 線性回歸線: y = 1.91526 * x + -2.29988
# 理論2階精度線通過參考點 (-2.99573, -8.04592)
f(x) = 1.91526 * x + -2.29988
g(x) = 2.0 * (x - -2.99573) + -8.04592

plot 'grid_convergence_data.dat' using 3:5 with linespoints pt 7 ps 1.5 lw 2 title sprintf('Computed (slope = %.2f)', 1.91526), \
     f(x) with lines lw 2 lc rgb 'red' title sprintf('Linear Fit (slope = %.2f)', 1.91526), \
     g(x) with lines lw 2 lc rgb 'green' dashtype 2 title '2nd Order Theory (slope = 2.0)'
