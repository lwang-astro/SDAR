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

mlst=`seq 0.05 0.05 1.0`
prefix='triple.m3.'
for m3 in $mlst
do
    echo '0 0 '$m12' '$m3' '$semi2' '$ecc2' '$inc2' '$roth2' '$rots2' '$ecca2 >$prefix$m3.orbit
    echo '1 0 '$m1'  '$m2' '$semi1' '$ecc1' '$inc1' '$roth1' '$rots1' '$ecca1 >>$prefix$m3.orbit
    echo 3 >$prefix$m3
    keplertree -n 2 $prefix$m3.orbit >>$prefix$m3
    echo '1 0 2 0 1' >>$prefix$m3

    echo $m3
    ./hermite -r 0.2 -t 100.0 --slowdown-ref '1e-6' ../input/$prefix$m3 >m$m3.log 2>m$m3.err
    ./hermite -r 0.2 -t 100.0 --slowdown-ref '5e-2' ../input/$prefix$m3 >m$m3.sd.log 2>m$m3.sd.err
done
