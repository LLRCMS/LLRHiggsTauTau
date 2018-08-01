#### check if AODSIM dataset has PU2017 

export INPUT="$1"
while IFS= read line
do
    export DATASET="${line}"

    if  [[ $DATASET == \#* ]]; 
    then
	if [[ $DATASET == \#\#* ]];
	then
	    export comment="${line}"
	    echo $comment
	
	else
	    continue
	fi
    elif [[ $DATASET = "="* ]];
    then
	continue
    elif [[ $DATASET == "" ]];
    then
	continue
    else
	echo "  ${DATASET}"


	export status=$(dasgoclient --query="dataset = ${DATASET} |grep dataset.status" ;)

	export parent=$(dasgoclient --query="parent dataset = ${DATASET}" ;)
	#echo "    $parent"
	if [[ $parent =~ "PU2017" && $status =~ "VALID" ]]
	then
	    echo -e "     \e[1;32mPU is corrected\e[0m"
	    
	    
	else 
	    echo -e "     \e[1;31mPU is not corrected\e[0m"
	fi
    fi	
done<"$INPUT"
