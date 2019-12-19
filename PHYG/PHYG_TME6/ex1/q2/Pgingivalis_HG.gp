set terminal png tiny size 800,800
set output "Pgingivalis_HG.png"
set size 1,1
set grid
unset key
set border 15
set tics scale 0
set xlabel "Pgingivalis_33277"
set ylabel "Pgingivalis_HG66"
set format "%.0f"
set mouse format "%.0f"
set mouse mouseformat "[%.0f, %.0f]"
if(GPVAL_VERSION < 5) { set mouse clipboardformat "[%.0f, %.0f]" }
set xrange [1:2354886]
set yrange [1:2441780]
set style line 1  lt 1 lw 3 pt 6 ps 1
set style line 2  lt 3 lw 3 pt 6 ps 1
set style line 3  lt 2 lw 3 pt 6 ps 1
plot \
 "Pgingivalis_HG.fplot" title "FWD" w lp ls 1, \
 "Pgingivalis_HG.rplot" title "REV" w lp ls 2
