m3=1.0
m12=0.5
m1=0.25
m2=0.25
semi2=1.0
ecc2=0.0
inc2=0.0
roth2=0.0
rots2=0.0
ecca2=0.0
semi1=0.1
ecc1=0.5
inc1=2.5
roth1=0.0
rots1=0.0
ecca1=1.5

prefix='triple.k.'
echo '0 0 '$m12' '$m3' '$semi2' '$ecc2' '$inc2' '$roth2' '$rots2' '$ecca2 >$prefix.orbit
echo '1 0 '$m1'  '$m2' '$semi1' '$ecc1' '$inc1' '$roth1' '$rots1' '$ecca1 >>$prefix.orbit
echo 3 >$prefix
keplertree -n 2 $prefix.orbit >>$prefix
echo '1 0 2 0 1' >>$prefix

klst='1 3 5 7 9 11 13 15 17 19 21 23 25 27 29'
for k in $klst
do
    echo $k
    ./hermite -r 0.2 -t 100.0 --slowdown-ref $k'e-3' ../input/$prefix >k$k.log 2>k$k.err
done
