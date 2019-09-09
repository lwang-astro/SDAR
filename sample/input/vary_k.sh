klst='1 3 5 7 9 11 13 15 17 19 21 23 25 27 29'
for k in $klst
do
    echo $k
    ./hermite -r 0.2 -t 100.0 --slowdown-ref $k'e-3' ../input/triple.stable >k$k.log 2>k$k.err
done
