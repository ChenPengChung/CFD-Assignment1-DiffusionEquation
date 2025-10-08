set terminal png enhanced size 800,600
set output '3.4orderCD-2order.png'
set title '3.4orderCD-2order : L1 Error vs Grid Spacing'
set xlabel 'log(dx)'
set ylabel 'log(L1 Error)'
set grid
set key left top
f3(x) = 3.92477 * x + 0.704053
plot '3.4orderCD-2order.dat' using 3:5 pt 13 ps 1.5  lc 'brown' title sprintf('4orderCD-2order (slope = %.2f)', 3.92477), \
     f3(x) with lines lw 2 lc 'brown' title sprintf('4orderCD-2order linear (slope = %.2f)', 3.92477), \

