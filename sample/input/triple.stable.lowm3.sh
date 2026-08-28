# try different versions of the AR code for a stable triple system
# need to compile AR sample code first

method_list='logh.ttl.sd logh.ttl logh.sd logh'
# see meaning of method in README.md
for suffix in $method_list
do
    ../AR/ar.$suffix -t 1.0e-3 -n 1000 triple.stable.lowm3 >triple.stable.lowm3.$suffix.log
done

