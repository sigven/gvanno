cp /Users/sigven/research/software/vcf2tsv/vcf2tsv.py gvanno/
cp /Users/sigven/research/docker/pcgr/src/pcgr/lib/annoutils.py gvanno/lib/
tar czvfh gvanno.tgz gvanno/
echo "Build the Docker Image"
TAG=`date "+%Y%m%d"`
docker build --no-cache -t sigven/gvanno:$TAG --rm=true .

