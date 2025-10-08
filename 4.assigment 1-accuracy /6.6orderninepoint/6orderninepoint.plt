set terminal png enhanced size 800,600
set output '6.ninepoint-6order.png'
set title '6.ninepoint-6order'
set xlabel 'log(dx)'
set ylabel 'log(L1 Error)'
set grid
set key left top
f(x) = 5.24533 * x + -6.10883
plot '6orderninepoint.dat' using 3:5  pt 6 ps 1.5 lc 'blue' title sprintf('6orderninepoint (slope = %.2f)', 5.24533), \
     f(x) with lines lw 2 lc 'blue' title sprintf('ninepoint-6orderlinear (slope = %.2f)', 5.24533), \
