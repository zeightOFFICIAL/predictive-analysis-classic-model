set terminal pngcairo enhanced font 'Arial,12'
set output 'gold_vs_PL=F_closing_price.png'
set title 'Gold vs Platinum (r = 0.367)'
set xlabel 'Gold Price (USD)'
set ylabel 'Platinum Price'
set grid
plot 'plot_PL=F_closing_price.dat' with points pt 7 ps 0.5 title ''
