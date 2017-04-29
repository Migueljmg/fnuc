set terminal pngcairo size 1200,600
set output "1.png"

set samples 1e6
set yrange [0:3500]

a=1095.7
b=1182.3
c=1074
d=945
e=844
f=104.5
g=32.7
h=106.6
i=70.2
j=17.33
k=8.56
m=11.36
n=142.8
o=2700
p=2443
q=80.2
r=3.57

c1=291.5
w1=6.6
c2=369.5
w2=6.4
c3=426.14
w3=11.21
c4=469.21
w4=14.79
c5=552.21
w5=9.4
lmin=694.53
c6=732.18
w6=13.45
c7=819.52
w7=13.69
lmax=878.7

l1=200
pp=680
qq=600
lmin=pp
lmax=880
z=920

f(x)=(x<l1)?0:(x<c1-w1)?a:(x<c1+w1)?a+(b-a)*(x-c1+w1)/(2*w1):(x<c2-w2)?b+(c-b)*(x-c1-w1)/(c2-w2-c1-w1):(x<c2+w2)?c+(d-c)*(x-c2+w2)/(2*w2):(x<c3-w3)?d+(e-d)*(x-c2-w2)/(c3-w3-c2-w2):(x<c3+w3)?e+(f-e)*(x-c3+w3)/(2*w3):(x<c4-w4)?f+(g-f)*(x-c3-w3)/(c4-w4-c3-w3):(x<c4+w4)?g+(h-g)*(x-c4+w4)/(2*w4):(x<c5-w5)?h+(i-h)*(x-c4-w4)/(c5-w5-c4-w4):(x<c5+w5)?i+(j-i)*(x-c5+w5)/(2*w5):(x<pp)?j+(k-j)*(x-c5-w5)/(pp-c5-w5):0

print d-(e-d)*(c2+w2)/(c3-w3-c2-w2)
print (e-d)/(c3-w3-c2-w2)

g(x)=(x<qq)?0:(x<lmin)?m:(x<c6-w6)?m+(n-m)*(x-lmin)/(c6-w6-lmin):(x<c6+w6)?n+(o-n)*(x-c6+w6)/(2*w6):(x<c7-w7)?o+(p-o)*(x-c6-w6)/(c7-w7-c6-w6):(x<c7+w7)?p+(q-p)*(x-c7+w7)/(2*w7):(x<lmax)?q+(r-q)*(x-c7-w7)/(lmax-c7-w7):(x<z+1)?r:0

#fit [l1:pp-1] f(x) "1.txt" u 1:2:3 via a,b,c,d,e,f,g,h,i,j,k,c1,w1,c2,w2,c3,w3,c4,w4,c5,w5
#fit [qq:z] g(x) "1.txt" u 1:2:3 via m,n,o,p,q,r,lmin,c6,w6,c7,w7,lmax

plot "1.txt" w yerrorbars title "Dados Experimentais" lc rgb "grey", [qq:z] g(x) lc rgb "black" lw 3 title "Ajuste Numerico", [l1:pp-1] f(x) lc rgb "black" lw 3 notitle
