export INPUT="$1"
while IFS= read line
do
    export DATASET="${line}"

    if  [[ $DATASET == \#* ]]; 
    then
	export comment="${line}"
	echo $comment
    elif [[ $DATASET = "="* ]];
    then
	continue
    elif [[ $DATASET == "" ]];
    then
	continue
    else
	echo "  ${DATASET}"
	if [[ ! $DATASET =~ "12Apr2018" &&  ! $DATASET =~ "31Mar2018" ]];
	then
	    echo "  !!! OUTDATED DATASET"	
	fi

	export status=$(dasgoclient --query="dataset = ${DATASET} |grep dataset.status" ;)

	if [[ $status =~ "VALID" ]]
	then
	    echo "     status = VALID"
	    export nevents=$(dasgoclient  --query="summary dataset = ${DATASET} |grep summary.nevents";)
	    echo "     nevents = ${nevents}"
	elif  [[ $status == "" ]];
	then
	    echo '     NOT FOUND'
	    
	elif [[ ! $status =~ "VALID" ]];
	then
	    echo '     status = PRODUCTION'
	    export nevents=$(dasgoclient  --query="summary dataset = ${DATASET} |grep summary.nevents";)
	    echo "     nevents = ${nevents}";

	fi
    fi	
done<"$INPUT"
