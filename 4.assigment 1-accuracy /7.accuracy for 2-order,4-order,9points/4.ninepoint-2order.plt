set terminal png enhanced size 800,600
set output '4.ninepoint-2order.png'
set title '4.ninepoint-2order : L1 Error vs Grid Spacing'
set xlabel 'log(dx)'
set ylabel 'log(L1 Error)'
set grid
set key left top
f4(x) = 3.90919 * x + 0.532055
plot '4.ninepoint-2order.dat' using 3:5  pt 8 ps 1.5 lc'red' title sprintf('ninepoint-2order (slope = %.2f)', 3.90919), \
     f4(x) with lines lw 2 lc 'red' title sprintf('ninepoint-2order linear (slope = %.2f)', 3.90919), \
     
