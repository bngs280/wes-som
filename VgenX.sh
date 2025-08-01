sudo docker run -it \
	-v /home/ubuntu/rgnex_som/rgenXSom_15022025/dockerdata/VEP_data_Cache:/usr/src/app/vepC \
	-v /home/ubuntu/rgnex_som/rgenXSom_15022025/dockerdata/vep_pluginDB:/usr/src/app/vepDB \
	-v /home/ubuntu/rgnex_som/rgenXSom_15022025/dockerdata/reference:/usr/src/app/ref17 \
	-v /home/ubuntu/rgnex_som/Dhit:/usr/src/app/input \
	-v /home/ubuntu/rgnex_som/Dhit/output/130325hg19:/usr/src/app/output \
	-v /home/ubuntu/rgnex_som/rgenXSom_15022025/nextflow.config:/usr/src/app/nextflow.config \
	-v /home/ubuntu/rgnex_som/rgenXSom_15022025/config.json:/usr/src/app/config.json \
	-v /home/ubuntu/rgnex_som/Dhit/output/130325hg19/work:/usr/src/app/work \
	rgenxtools:latest
