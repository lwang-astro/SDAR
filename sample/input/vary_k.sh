# check how slowdown factor affect the evolution of a triple system by using hermite code
# need to compile hermite and keplertree first

# mass 
m3=1.0
m12=0.5 # m1+m2
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
semi1=0.1  # semi-major axis             
ecc1=0.5   # eccentricity                
inc1=2.5   # inclination                 
roth1=0.0  # rotation angle in x-y plane 
rots1=0.0  # rotation angle in rest frame
ecca1=1.5  # eccentricty anomaly         

# use keplertree to generate particle data from orbital data
prefix='triple.k'
echo '0 0 '$m12' '$m3' '$semi2' '$ecc2' '$inc2' '$roth2' '$rots2' '$ecca2 >$prefix.orbit
echo '1 0 '$m1'  '$m2' '$semi1' '$ecc1' '$inc1' '$roth1' '$rots1' '$ecca1 >>$prefix.orbit
../Kepler/keplertree -n 2 $prefix.orbit > $prefix.particle

# generate input for hermite
echo 3 >$prefix
awk '{print $LINE,0.0}' $prefix.particle >>$prefix
echo '1 0 2 0 1' >>$prefix

# slowdown factor parameter 
klst='1 3 5 7 9 11 13 15 17 19 21 23 25 27 29'
for k in $klst
do
    echo $k
    ../Hermite/hermite -r 0.2 -t 100.0 --slowdown-ref $k'e-3' $prefix >k.$k.log 2>k.$k.err
done
