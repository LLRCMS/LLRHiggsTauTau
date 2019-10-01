#Ensure you source crab and setup voms certificate
for filename in crab3/crab_Oct19_*/; do
  echo $filename
  crab resubmit $filename
done
