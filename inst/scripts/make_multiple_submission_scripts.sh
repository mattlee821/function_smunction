# create multiple submission scripts from a single master script

## make a master script = master_qsub

## under the bash variable section add a specified variabe = "VAR1=""
ls -1 > filenames # create a file with the names of all files within a directory

## create multiple .sh qsub scripts from a single file with names in based on a master_sub script in which line (NR == ??) holds the variable name
cat filenames | while read i; do echo ${i}; awk '{ if (NR == 17) print "VAR1='${i}'"; else print $0}' master_qsub > ${i}.sh; done 

## you can create a single submission script which will run all of the jobs in one go 
ls |grep FILE |sed 's/\ / /g' |awk '{print "qsub "$i" -e log -o log"}' > run.sh # replace FILE with the common character for your scripts you made in the previosu step
sh run.sh
