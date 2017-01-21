tar czvfh gvanno.tgz gvanno/
echo "Build the Docker Image"
TAG=`date "+%Y%m%d"`
docker build -t sigven/gvanno:$TAG --rm=true .

