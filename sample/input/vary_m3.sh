# check how third body mass affects the evolution of a triple system by using hermite code
# need to compile hermite and keplertree first

#mass
m12=0.5 #m1+m2
m1=0.25
m2=0.25

# outer orbit (G=1)
semi2=1.0 # semi-major axis             
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

# m3 using a sequence 
mlst=`seq 0.05 0.05 1.0`

prefix='triple.m3'
for m3 in $mlst
do
    # use keplertree to generate particle data from orbital data
    echo '0 0 '$m12' '$m3' '$semi2' '$ecc2' '$inc2' '$roth2' '$rots2' '$ecca2 >$prefix.$m3.orbit
    echo '1 0 '$m1'  '$m2' '$semi1' '$ecc1' '$inc1' '$roth1' '$rots1' '$ecca1 >>$prefix.$m3.orbit
    ../Kepler/keplertree -n 2 $prefix.$m3.orbit >>$prefix.$m3.particle

    # generate input
    echo 3 >$prefix.$m3
    awk '{print $LINE,0.0}' $prefix.$m3.particle >>$prefix.$m3
    echo '1 0 2 0 1' >>$prefix.$m3

    echo $m3
    ../Hermite/hermite -r 0.2 -t 100.0 --slowdown-ref '1e-6' $prefix.$m3 >m3.$m3.log 2>m3.$m3.err
    ../Hermite/hermite -r 0.2 -t 100.0 --slowdown-ref '5e-2' $prefix.$m3 >m3.$m3.sd.log 2>m3.$m3.sd.err
done
