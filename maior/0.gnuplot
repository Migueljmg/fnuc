set terminal pngcairo size 1200,600
set output "0.png"

set samples 1e6
set yrange [0:3500]

l1=220
a=1349.1
b=1408
c=1375.5
d=1085.2
e=77.9
f=37.9
g=112.9
h=81.6
i=18.3
j=8.24
m=11.2
n=183
o=2945
p=2663
q=117.6
r=7.3

c1=310
w1=6
c2=370.6
w2=23.3
c3=444.73
w3=13.53
c4=488.07
w4=10.22
c5=553.71
w5=8.6
lmin=718.54
c6=751.98
w6=11.6
c7=824.25
w7=10.78
ldir=871.7

pp=680
fin=900

z=fin-1
f(x)=(x<l1)?0:(x<c1-w1)?a:(x<c1+w1)?a+(b-a)*(x-c1+w1)/(2*w1):(x<c2-w2)?b+(c-b)*(x-c1-w1)/(c2-w2-c1-w1):(x<c2+w2)?c+(d-c)*(x-c2+w2)/(2*w2):(x<c3-w3)?d:(x<c3+w3)?d+(e-d)*(x-c3+w3)/(2*w3):(x<c4-w4)?e+(f-e)*(x-c3-w3)/(c4-w4-c3-w3):(x<c4+w4)?f+(g-f)*(x-c4+w4)/(2*w4):(x<c5-w5)?g+(h-g)*(x-c4-w4)/(c5-w5-c4-w4):(x<c5+w5)?h+(i-h)*(x-c5+w5)/(2*w5):(x<pp)?i+(j-i)*(x-c5-w5)/(pp-c5-w5):0

g(x)=(x<pp)?0:(x<lmin)?m:(x<c6-w6)?m+(n-m)*(x-lmin)/(c6-w6-lmin):(x<c6+w6)?n+(o-n)*(x-c6+w6)/(2*w6):(x<c7-w7)?o+(p-o)*(x-c6-w6)/(c7-w7-c6-w6):(x<c7+w7)?p+(q-p)*(x-c7+w7)/(2*w7):(x<ldir)?q+(r-q)*(x-c7-w7)/(ldir-c7-w7):(x<fin)?r:0

#fit [l1:z] f(x) "0.txt" u 1:2:3 via a,b,c,d,e,f,g,h,i,j,c4,w4,c5,w5
#fit [pp:z] g(x) "0.txt" u 1:2:3 via m,n,o,p,q,r,lmin,c6,w6,c7,w7,ldir

plot "0.txt" w yerrorbars title "Dados Experimentais" lc rgb "grey", [pp:z] g(x) lc rgb "black" lw 3 title "Ajuste Numerico", [l1:pp] f(x) lc rgb "black" lw 3 notitle
