# Set the output format to EPS with enhanced font support
set terminal postscript eps enhanced color font 'Helvetica,18'

# Set the output filename
set output 'final_new_currents_plot.eps'

# Set line and point properties
set style line 1 lc rgb 'black' lw 3             # Solid line for currents.txt
set style line 2 lc rgb 'black' lw 3             # Solid line for currentsMoM.txt
set style line 3 lc rgb 'black' pt 2 ps 2 # 'x' marker for currents.txt
set style line 4 lc rgb 'black' pt 1 ps 3 # '+' marker for currentsMoM.txt
set style line 5 lc rgb 'black' lw 3 pt 2 ps 2 
set style line 6 lc rgb 'black' lw 3 pt 1 ps 3

# Set labels, font sizes, and other plot properties
set xlabel "z(m)" font 'Helvetica,16'
set ylabel "Currents (A)" font 'Helvetica,16'
set title "Conventional MoM vs Method of Forces" font 'Helvetica,18'
set grid
set key font 'Helvetica,18'
set xtics font 'Helvetica,18'
set ytics font 'Helvetica,18'

# Plot the data
plot "currentsMoF.txt" using 1:2 with lines linestyle 1 notitle , \
     "currentsMoM.txt" using 1:2 with lines linestyle 2 notitle  ,\
     keyentry with linespoints ls 5 title "MoF",\
     keyentry with linespoints ls 6 title "MoM",\
     "currentsMoF.txt" using 1:2 every 4::1::41 with points  linestyle 3 notitle , \
     "currentsMoM.txt" using 1:2 every 4::3::41 with points linestyle 4 notitle 

# End the script
quit

