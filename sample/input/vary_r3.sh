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

alst=`seq 1.0 0.2 5.0`
prefix='triple.a2.'
for a2 in $alst
do
    echo '0 0 '$m12' '$m3' '$a2' '$ecc2' '$inc2' '$roth2' '$rots2' '$ecca2 >$prefix$a2.orbit
    echo '1 0 '$m1'  '$m2' '$semi1' '$ecc1' '$inc1' '$roth1' '$rots1' '$ecca1 >>$prefix$a2.orbit
    echo 3 >$prefix$a2
    keplertree -n 2 $prefix$a2.orbit >>$prefix$a2
    echo '1 0 2 0 1' >>$prefix$a2

    echo $a2
    ./hermite -r 0.2 -t 100.0 --slowdown-ref '1e-9' ../input/$prefix$a2 >a$a2.log 2>a$a2.err
    ./hermite -r 0.2 -t 100.0 --slowdown-ref '5e-2' ../input/$prefix$a2 >a$a2.sd.log 2>a$a2.sd.err
done
