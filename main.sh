sudo docker run -it \
	-v /home/ubuntu/rgnex_som/wes_som_data/dockerdata/VEP_data_Cache:/usr/src/app/vepC \
	-v /home/ubuntu/rgnex_som/wes_som_data/dockerdata/vep_pluginDB:/usr/src/app/vepDB \
	-v /home/ubuntu/rgnex_som/wes_som_data/dockerdata/reference:/usr/src/app/ref17 \
	-v /home/ubuntu/rgnex_som/wes_som_out:/usr/src/app/input \
	-v /home/ubuntu/rgnex_som/wes_som_out/output/130325hg19:/usr/src/app/output \
	-v /home/ubuntu/rgnex_som/wes_som_data/nextflow.config:/usr/src/app/nextflow.config \
	-v /home/ubuntu/rgnex_som/wes_som_data/config.json:/usr/src/app/config.json \
	-v /home/ubuntu/rgnex_som/Dhit/output/130325hg19/work:/usr/src/app/work \
	wes_som:latest
