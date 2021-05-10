
gcc -fopenmp -O3 -o driver driver.c

declare -a threads=( 1 2 5 10 20 40 60 )

counter2=256

while [ $counter2 -le 8192 ]
do
    for i in "${threads[@]}"
    do
        export OMP_NUM_THREADS=$i
        ex1=1
        ex2=1
        if [ $i -gt 1 ]
        then
            ex1=0
            ex2=0
        fi
        echo -n "$counter2, $i, "
        echo -n "$counter2, $i, " >>result.txt
        ./driver $counter2 $ex1 $ex2 1 1 1
        ./driver $counter2 $ex1 $ex2 1 1 1 >>result.txt
    done
    counter2=$(( $counter2 * 2 ))
done
