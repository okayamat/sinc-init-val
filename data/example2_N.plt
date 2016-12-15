set encoding iso_8859_1
set term postscript eps enhanced "Times-Roman" 20
set output "example2_N.eps"
set size 0.82
set logscale y
set key spacing 3
set xrange [0:140]
set yrange [1e-16:1]
set xlabel "{/Times-Italic=24 N}"
set ylabel "{/Times-Roman=24 maximum error}"
plot "SE_nyst_example2.dat" using 1:3 w lp title "SE-Sinc-Nystr{\366}m" lt 3 pt 2 ps 2, "SE_coll_example2.dat" using 1:3 w lp title "SE-Sinc-collocation" lt 3 pt 8 ps 2, "DE_nyst_example2.dat" using 1:3 w lp title "DE-Sinc-Nystr{\366}m" lt 4 pt 1 ps 2, "DE_coll_example2.dat" using 1:3 w lp title "DE-Sinc-collocation" lt 4 pt 6 ps 2
