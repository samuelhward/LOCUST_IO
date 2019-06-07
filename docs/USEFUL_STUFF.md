## Useful Stuff:

* [Fortran format guide](https://pages.mtu.edu/~shene/COURSES/cs201/NOTES/format.html)
* [GEQDSK format](http://nstx.pppl.gov/nstx/Software/Applications/a-g-file-variables.txt)
* get_fbm guide - ftp://ftp.pppl.gov/pub/transp/codata/doc/fbm.doc
* [TRANSP help](https://w3.pppl.gov/~pshare/help/transp.htm)

## Handy Functions & Aliases:

```bash
    function lwatch(){ #first arg is tokhead, second is number of preview lines, third is the run_data file number
        watch tail -n "$2" "/tmp/${USER}/locust.$1/OutputFiles/*$3"
   }

    function locust_mkdir(){
        mkdir /tmp/${USER}/
        for tokhead in "$@"
        do
            mkdir /tmp/${USER}/locust.$tokhead/
            mkdir /tmp/${USER}/locust.$tokhead/InputFiles
            mkdir /tmp/${USER}/locust.$tokhead/OutputFiles
            mkdir /tmp/${USER}/locust.$tokhead/CacheFiles
        done
    }

    function locust_pull(){
        for version in "$@"
        do
            (cd ~/locusts/locust_$version/locust/ && make clean && git stash && git pull)
        done
    }

    function flag_grab(){
        FLAGS=$(grep '\-D' "$1" | awk '{ print $1 }' | grep '\-D' | paste -sd " ")

        if [ -z "$2" ]
            then
            echo "make FLAGS='${FLAGS}'" #no second arg given, print FLAGS to screen
        else
            echo "make FLAGS='${FLAGS}'" > $2 #if second arg given, echo to that filename
        fi
    }
```