set terminal pngcairo size 1200,600
set output "reta.png"

set samples 1e6

f(x)=a+(b*10**(-3))*x
g(y)=(y-a)/(b*10**(-3))

fit g(x) "reta.txt" u 2:1:3 via a,b

plot "reta.txt" u 1:2:3:(1e-8) w xyerrorbars title "Pontos de Calibracao" lc rgb "grey", f(x) lc rgb "black" title "Reta de Calibracao"
