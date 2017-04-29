set terminal pngcairo size 1200,600
set output "calib.png"

set samples 1e6
set yrange [0:1600]

l1=240
l2=290
l3=314
l4=493
l5=524
l6=730
l7=763
l8=854
l9=870
l10=885

z=l10-1

c1=(l3+l2)/2
w1=(l3-l2)/2
c2=(l5+l4)/2
w2=(l5-l4)/2
c3=(l7+l6)/2
w3=(l7-l6)/2
c4=(l9+l8)/2
w4=(l9-l8)/2

f(x)=(x<l1)?0:(x<l2)?a:(x<l3)?a+b*(x-l2):0
f(x)=(x<l1)?0:(x<c1-w1)?a:(x<c1+w1)?a+(b-a)*(x-c1+w1)/(2*w1):(x<c2-w2)?b:(x<c2+w2)?b+(c-b)*(x-c2+w2)/(2*w2):(x<c3-w3)?c:(x<c3+w3)?c+(d-c)*(x-c3+w3)/(2*w3):(x<c4-w4)?d:(x<c4+w4)?d+(e-d)*(x-c4+w4)/(2*w4):(x<l10)?e:0

fit [l1:z] f(x) "calib.txt" u 1:2:3 via a,b,c,d,e,c1,w1,c2,w2,c3,w3,c4,w4

plot "calib.txt" w yerrorbars title "Dados Experimentais" lc rgb "grey", [l1:z] f(x) lc rgb "black" lw 3 title "Ajuste Numerico"
