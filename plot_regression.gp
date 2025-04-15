set terminal pngcairo enhanced font 'Arial,12'
set output 'silver_gold_regression.png'
set title 'Silver vs Gold Regression (R² = 0.742845)'
set xlabel 'Gold Price (USD)'
set ylabel 'Silver Price (USD)'
set grid
plot 'silver_gold_regression.dat' using 1:2 with points pt 7 ps 0.5 lc rgb 'blue' title 'Actual', \
     'silver_gold_regression.dat' using 1:3 with lines lw 2 lc rgb 'red' title 'Fitted'
