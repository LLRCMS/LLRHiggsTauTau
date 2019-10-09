#Ensure you source crab and setup voms certificate
for filename in crab3_Oct19/crab_test_Oct19*/; do
  echo $filename
  crab resubmit -d $filename
done
