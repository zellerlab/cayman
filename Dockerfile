FROM ubuntu:22.04

LABEL maintainer="cschu1981@gmail.com"
LABEL version="0.10.1"
LABEL description="cayman - profiling carbohydrate active enzymes in metagenomic/transcriptomic wgs samples"


ARG DEBIAN_FRONTEND=noninteractive

RUN apt update
RUN apt upgrade -y

RUN apt install -y wget python3-pip git dirmngr gnupg ca-certificates build-essential libssl-dev libcurl4-gnutls-dev libxml2-dev libfontconfig1-dev libharfbuzz-dev libfribidi-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev bwa
RUN apt clean

RUN mkdir -p /opt/software && \
	cd /opt/software && \
	git clone https://github.com/zellerlab/cayman && \
	cd cayman && \
	pip install .
  
CMD ["cayman"]
