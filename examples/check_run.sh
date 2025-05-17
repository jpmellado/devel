#/bin/bash

if [ $# -eq 0 ]; then
    echo "Usage: $0 bin-directory [mpirun-command bin-directory-preproc]"
    exit 1
fi

BINPATH=$1
MPIRUN=$2
if [ $# -eq 3 ]; then
    BINPATH_PP=$3
    MPIRUN_PP=""
else
    BINPATH_PP=$BINPATH
    MPIRUN_PP=$MPIRUN
fi

if [ -e "dns.out" ]; then
    echo -e "\033[1;31mFailed \033[0m[dns.out exists]."
    exit 1
else

    if [ -e "tlab.ini" ]; then

        #PreProcessing
        $MPIRUN_PP $BINPATH_PP/inigrid.x
        if [ $? = 0 ]; then
            $MPIRUN_PP $BINPATH_PP/inirand.x
            if [ $? = 0 ]; then
                $MPIRUN_PP $BINPATH_PP/iniscal.x
                if [ $? = 0 ]; then
                    $MPIRUN_PP $BINPATH_PP/iniflow.x
                    # if [ $? = 0 ]; then
                    #     $MPIRUN $BINPATH/inipart.x

                    #Simulation
                    if [ $? = 0 ]; then
                        LIST=$(ls *.ics*)
                        for FILE in $LIST; do mv $FILE ${FILE/ics/0}; done

                        $MPIRUN $BINPATH/dns.x
                        if [[ $? = 0 && ! -e "tlab.err" ]]; then
                            diff dns.out dns.out.ref >/dev/null 2>&1
                            if [ $? = 0 ]; then
                                grep -i " nan " avg* >/dev/null 2>&1
                                if [ $? = 0 ]; then
                                    echo -e "\033[1;31mFailed \033[0m[NaN in averages]."
                                    exit 10
                                else
                                    echo -e "\033[1;32mPassed\033[0m."
                                fi
                            else
                                echo -e "\033[1;31mFailed \033[0m[dns.out]."
                                exit 9
                            fi
                        else
                            echo -e "\033[1;31mFailed \033[0m[dns]."
                            exit 8
                        fi
                        # else
                        # echo -e "\033[1;31mFailed \033[0m[inipart]."; exit 7
                        # fi
                    else
                        echo -e "\033[1;31mFailed \033[0m[iniflow]."
                        exit 6
                    fi
                else
                    echo -e "\033[1;31mFailed \033[0m[[iniscal]."
                    exit 5
                fi
            else
                echo -e "\033[1;31mFailed \033[0m[[inirand]."
                exit 4
            fi
        else
            echo -e "\033[1;31mFailed \033[0m[[inigrid]."
            exit 3
        fi

    else
        echo -e "\033[1;31mFailed \033[0m[[tlab.ini]."
        exit 2
    fi

fi
