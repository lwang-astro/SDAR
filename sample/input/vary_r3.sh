# check how outer semi-major axis affects the evolution of a triple system by using hermite code
# need to compile hermite and keplertree first

#mass
m3=1.0
m12=0.5 #m1+m2
m1=0.25
m2=0.25

# outer orbit (G=1)
ecc2=0.0  # eccentricity                
inc2=0.0  # inclination                 
roth2=0.0 # rotation angle in x-y plane 
rots2=0.0 # rotation angle in rest frame
ecca2=0.0 # eccentricty anomaly         

# inner orbit (G=1)
semi1=0.1 # semi-major axis             
ecc1=0.5  # eccentricity                
inc1=2.5  # inclination                 
roth1=0.0 # rotation angle in x-y plane 
rots1=0.0 # rotation angle in rest frame
ecca1=1.5 # eccentricty anomaly         

semi2lst=`seq 1.0 0.2 5.0`
prefix='triple.a2'
for semi2 in $semi2lst
do
    # use keplertree to generate particle data from orbital data
    echo '0 0 '$m12' '$m3' '$semi2' '$ecc2' '$inc2' '$roth2' '$rots2' '$ecca2 >$prefix.$semi2.orbit
    echo '1 0 '$m1'  '$m2' '$semi1' '$ecc1' '$inc1' '$roth1' '$rots1' '$ecca1 >>$prefix.$semi2.orbit
    ../Kepler/keplertree -n 2 $prefix.$semi2.orbit >>$prefix.$semi2.particle

    # generate input
    echo 3 >$prefix.$semi2
    awk '{print $LINE,0.0}' $prefix.$semi2.particle >>$prefix.$semi2
    echo '1 0 2 0 1' >>$prefix.$semi2

    echo $semi2
    ../Hermite/hermite -r 0.2 -t 100.0 --slowdown-ref '1e-9' $prefix.$semi2 >semi2.$semi2.log 2>semi2.$semi2.err
    ../Hermite/hermite -r 0.2 -t 100.0 --slowdown-ref '5e-2' $prefix.$semi2 >semi2.$semi2.sd.log 2>semi2.$semi2.sd.err
done
