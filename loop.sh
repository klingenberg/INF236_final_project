
gcc -fopenmp -O3 -o driver driver.c

declare -a threads=( 1 2 5 10 20 40 60 )

counter2=256

while [ $counter2 -le 4096 ]
do
    counter=2
    while [ $counter -le 256 ]
    do
        for i in "${threads[@]}"
        do
            export OMP_SET_NUM_THREADS=$i
            ex1=1
            ex2=1
            if [ $i -gt 1 ]
            then
                ex1=0
                ex2=0
            fi
            echo $counter2 $counter $i $ex1 $ex2 1 1 1
            ./driver $counter2 $counter $ex1 $ex2 1 1 1
        done
        counter=$(( $counter * 2 ))
    done
    counter2=$(( $counter2 * 2 ))
done
