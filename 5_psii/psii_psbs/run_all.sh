function run_all () {
    set -e
    sim=(2 3 4 5 6 7 8)
    for n_sim in "${sim[@]}"; do
        bash simulation.sh ${n_sim} > /dev/null 2>&1 &
    done

}
run_all