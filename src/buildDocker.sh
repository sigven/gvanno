#cp /Users/sigven/research/docker/pcgr/src/loftee_1.0.3.tgz .
#cp /Users/sigven/research/software/vcf2tsv/vcf2tsv.py gvanno/
tar czvfh gvanno.tgz gvanno/
echo "Build the Docker Image"
TAG=`date "+%Y%m%d"`
docker build --no-cache -t sigven/gvanno:$TAG --rm=true .

